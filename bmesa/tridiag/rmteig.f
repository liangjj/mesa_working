*deck rmteig.f
c***begin prologue     rmteig
c***date written       950704   (yymmdd)
c***revision date               (yymmdd)
c***keywords           r-matrix, diagonalization
c***author             schneider, barry (nsf)
c***source             mesa
c***purpose            calculate eigenvalues and eigenfunctions of 1-d
c***                   schroedinger equation by finite difference.         
c***   
c***                   
c***references         
c
c***routines called 
c***                   
      program rmteig
      implicit integer(a-z)
      parameter ( nrmax=100)
      common ia(1)
      real*8 z(1), rleft, rright, stp, dele, fpkey, err
      dimension npts(nrmax,2), nwts(nrmax,2), pt(2), wt(2)
      dimension points(2), rleft(nrmax), rright(nrmax), stp(2)
      character*4096 ops
      character*128 fillam
      character*16 pottyp, bc
      character*24 cpass
      character*1600 card
      character*80  chrkey, title, drct, typefd
      logical logkey, prntv, prntg, prntsl, prnth
      logical posinp, evecs, extrap, prntvc
      equivalence (ia(1),z(1))
      common /memory/ ioff
      common/io/inp,iout
c
      call drum
      call iosys ('read character options from rwf',-1,0,0,ops)
      call iosys ('read character "linear algebraic filename" from rwf',
     1                                                  -1,0,0,fillam)
      call iosys ('open lamdat as unknown',0,0,0,fillam)
      write(iout,*)
      title='r-matrix functions via finite difference'
      write(iout,*) title
      write(iout,*)
      prntv=logkey(ops,'print=m6238=potential',.false.,' ')
      prntg=logkey(ops,'print=m6238=grid',.false.,' ')
      prntsl=logkey(ops,'print=m6238=solution',.false.,' ')
      prnth=logkey(ops,'print=m6238=hamiltonian',.false.,' ')
      prntvc=logkey(ops,'print=m6238=vectors',.false.,' ')
      bc=chrkey(ops,'boundary-condition','zero-function',' ')
      l=intkey(ops,'angular-momentum',0,' ')
      pottyp=chrkey(ops,'potential-type','exponential',' ')
      evecs=logkey(ops,'eigenvectors',.false.,' ')
      order=intkey(ops,'order',3,' ')
      iter=intkey(ops,'maximum-iterations',200,' ')
      dele=fpkey(ops,'fractional-bound',1.d-01,' ')
      err=fpkey(ops,'relative error',1.d-10,' ')
      typefd=chrkey(ops,'difference-approximation',
     1                  'standard-finite-difference',' ')
      drct=chrkey(ops,'diagonalize','all',' ')
      if(drct.ne.'all') then
         evecs=.true.
      endif         
      extrap=logkey(ops,'extrapolate',.false.,' ')
c     basic coarse grid information and region definition
      if (posinp('$finite-difference',cpass) ) then
          call cardin(card)
          call dtafd(card,rleft,rright,nreg,npts,nwts,pt,wt,
     1               points,nrmax)
      endif
      nn=points(1)-2
      if(drct.ne.'all') then
         nroots=intkey(ops,'number-of-roots',5,' ')
      endif
      nroots=min(nn,nroots)
      bw=1
      if (order.eq.5) then
          bw=2
      endif
      if (pt(1).eq.0) then
          call lnkerr('error in point allocation')
      endif
      write(iout,*) title                    
      write(iout,1)  points(1), pottyp, l, bc, nroots
      write(iout,2) rleft(1), rright(nreg)
      dim=2
      if(.not.extrap) then
         write(iout,3)
         dim=1
         pt(2)=pt(1)
         wt(2)=wt(1)
         points(2)=points(1)
      endif
c     generate basic grid information
      ioff=1
      do 10 i=1,2
c        storage for all grids, potentials, effective potentials
c        and inhomogeneities are outlayed
         rpt=ioff
         wts=rpt+pt(2)
         sumwts=wts+wt(2) 
         runiqe=sumwts+pt(2) 
         v=runiqe+pt(2)
         band=v+pt(2)
         work=band+pt(2)*order
         dum=work+max(pt(2)*order,pt(2)*pt(2))
         ind=wpadti(dum+6*pt(2))
         eig=iadtwp(ind+pt(2)+pt(2))
         evec=eig+pt(1)+pt(2)
         s=evec+pt(1)*pt(1)+pt(2)*pt(2)
         wdsusd=s+pt(2)*pt(2)
         if (i.eq.1) then
             need=wpadti(wdsusd)
             call getscm(need,z,maxcor,'finite-difference',0)
         endif
 10   continue
      write(iout,*) 
      title='generating grid and potential for finite'
      len=length(title)
      title=title(1:len)//' difference region'
      write(iout,*) title
      ei=eig
      ev=evec
      do 20 i=1,dim
         call mkgrd(z(rpt),z(runiqe),z(wts),z(sumwts),rleft,rright,
     1              stp(i),nreg,npts(1,i),nrmax,prntg)
         call potntl(z(v),z(runiqe),points(i),pottyp,l,prntv)
         if (typefd.eq.'standard-finite-difference'
     1                         .or.
     2       typefd.eq.'standard-numerov') then
             call fdiff(z(band),z(v),stp(i),points(i),order,bw,typefd)
         else
             call gfdiff(z(band),z(runiqe),z(v),points(i),order,
     1                   bw,typefd)
         endif
         call bndcnd(z(band),points(i),order,bw,bc)
         call hdiag(z(band),z(ei),z(ev),z(work),z(dum),z(s),ia(ind),
     1              points(i),order,bw,evecs,nroots,prntvc)
         if(evecs) then
            call addend(z(ev),nroots,points(i),bc)
            call orth(z(s),z(ev),z(wts),z(runiqe),z(work),
     1                nroots,points(i))
         endif
         ei=ei+pt(i)
         ev=ev+pt(i)*pt(i)
 20   continue   
c     use richardson extrapolation technique
      if(extrap) then
         call richdn(z(eig),z(eig+pt(1)),z(evec),z(evec+pt(1)*pt(1)),
     1               stp,points(1),points(2),nroots)
      else
c
c     find better eigenvalues and eigenfunctions using these as guesses
c     and numerov method.
c
         call better(z(band),z(runiqe),z(v),z(eig),z(evec),
     1               z(s),z(wts),ia(ind),z(work),stp(1),nroots,
     2               points(1),order,bw,iter,dele,err,bc,
     3               typefd,prntsl)
      endif
      call chainx(0)
c
 1    format(/,1x,'number of finite difference points = ',i4,/,1x,
     1            'potential type                     = ',a16,/,1x,
     2            'angular momentum                   = ',i2,/1x,
     3            'type boundary condition            = ',a16,/,1x,
     4            'number of roots                    = ',i3)
 2    format(/,1x,'finite difference region begins at = ',e15.8,/,1x,
     1            'finite difference region ends at   = ',e15.8)
 3    format(/,1x,'do not perform richardson extrapolation')
      stop
      end

