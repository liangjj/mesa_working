*deck hcoll.f 
c***begin prologue     hcoll
c***date written       020303   (yymmdd)
c***revision date               (yymmdd)
c***keywords           integral equation
c***                   
c***author             schneider, b. i.(nsf)
c***source             hcoll
c***purpose            driver electron + h atom collisions
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       hcoll
      program hcoll
c
      implicit integer (a-z)
      parameter ( nenmax=100)
      character*4096 ops
      character*2 itoc
      character*800 card
      logical dollar, logkey
      logical prn, exch, opt 
      real*8 energy, fpkey, rmax, stp, g, psi_int
      dimension prn(10), energy(nenmax), wds(10), psi_int(500)
      common/io/inp, iout      
      pointer (pg,g(1)), (pg,ig(1))
c
      call drum
      write(iout,*)
      call iosys ('read character options from rwf',-1,0,0,ops)
      opt=logkey(ops,'optical-potential',.false.,' ')
      exch=logkey(ops,'exchange',.false.,' ')
      nen=intkey(ops,'number-of-energies',1,' ')
      call fparr(ops,'energies',energy,nen,' ')
      rmax=fpkey(ops,'maximum-r-value',20.d0,' ')
      npts=intkey(ops,'number-of-mesh-points',100,' ')
      e_1s=fpkey(ops,'target-energy',-.5d0,' ')
      stp=rmax/npts
      npts=npts+1
      nsol=2+nopt

c     Keep all the solutions to the equations contiguous
c

      psi_00=1
      psi_0=psi_00+npts
      words=psi_0+npts
      if(exch) then
         psi_1=words
         words=psi_1+npts
      endif
      if(opt) then
         psi_lam=words
         words=psi_lam+nopt*npts
      endif
      r=words
      r_k=r+npts
      i_k=r_k+npts
      v_loc=i_k+npts
      k_hat=v_loc+npts
      words=k_hat+npts+2*npts
      if(exch) then
         psi_1s=words
         vpsi_1s=psi_1s+npts
         u_hat=vpsi_1s+npts
         words=u_hat+2*npts
      endif
      if(opt) then
         v_lam=words
         e_lam=v_lam+nopt*npts
         words=e_lam+nopt
      endif
      if(exch.or.opt) then
         driver=words
         cmat=driver+npts
         c_0=cmat+nsol*nsol
         ipvt=wpadti(c_0+nsol)
         words=iadtwp(ipvt+nsol)
      endif
      words=wpadti(words)
      call memory(words,pg,wds(1),'main',0)
      call grid(g(r),stp,npts)
      call local(g(v_loc),g(r),npts)
      if(exch) then
         call f1s(g(psi_1s),g(vpsi_1s),g(r),npts)
      endif
      if(opt) then
         call rdlam(g(v_lam),g(e_lam),nopt,npts)
      endif
      do 10 i=1,nen
         call green(g(r_k),g(i_k),g(r),energy(i),npts)

c
c        solve for $\Psi_{00}$
c
         if(exch.or.opt) then
            call exnonit(g(psi_00),g(r),g(psi_1s),g(r_k),g(r_k),g(i_k),
     1                   g(v_loc),g(k_hat),g(u_hat),g(vpsi_1s),stp,
     2                   psi_int(1),,exch,npts)
            call phys(g(psi_0),g(psi_00),g(psi_00),
     1                psi_int(1),psi_int(1),npts)
c        solve with $\psi_{1s}$ as the inhomogeneity

            call drv(g(driver),g(psi_1s),g(r_k),g(i_k),stp,npts)
            call exnonit(g(psi_1),g(r),g(psi_1s),g(driver),g(r_k),
     1                   g(i_k),g(v_loc),g(k_hat),g(u_hat),
     2                   g(vpsi_1s),stp,psi_int(2),npts)
            call phys(g(psi_1),g(psi_00),g(psi_1),
     1                psi_int(1),psi_int(2),npts)
            if(opt) then
               int=3
               vlam=v_lam
               psilam=psi_lam
               do 20 j=1,nopt
                  call drv(g(driver),g(vlam),g(r_k),g(i_k),stp,npts)
                  call exnonit(g(psilam),g(r),g(psi_1s),g(driver),
     1                         g(r_k),g(i_k),g(v_loc),g(k_hat),
     2                         g(u_hat),g(vpsi_1s),stp,psi_int(int),
     3                         npts)
                  call phys(g(psilam),g(psi_00),g(psilam),
     1                      psi_int(1),psi_int(int),npts)
                  int=int+1
                  vlam=vlam+npts
                  psilam=psilam+npts
 20            continue
            endif

c        Now that we have all the elementary solutions, we can solve for the 
c        physical solution

            call fulsol(g(psi_0),g(psi_1s),g(vpsi_1s),g(v_lam),g(e_lam),
     1                  g(cmat),g(c_0),ig(ipvt),energy(i),e_1s,
     2                  stp,nopt,npts)
         else
            call nonit(g(psi_00),g(r),g(r_k),g(r_k),g(i_k),
     1                 g(v_loc),g(k_hat),stp,psi_int(1),npts)
            call phys(g(psi_0),g(psi_00),g(psi_00),
     1                psi_int(1),psi_int(1),npts)
         endif
 10   continue   

      call memory(-wds(1),pg,idum,'main',idum)
      call chainx(0)               
      stop
      end

















