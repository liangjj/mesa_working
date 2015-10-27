*deck psit.f 
c***begin prologue     psit
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, orthogonal polynomial
c***                   bose-einstein
c***author             schneider, b. i.(nsf)
c***source             tprop
c***purpose            time dependent schroedinger equation 
c***                   using time based dvr.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       psit
      subroutine psit(q,qm,wtm,pm,ham,v,trn,eig,
     1                pbasis,parray,spdim,dim,nptt,nmaxx,ug,mgr,spac,
     2                tim,timh0,i0stat,sptimv,state,coord,itsolv,
     3                cnverg,thresh,eps,precon,nblck,card,prnt,hdiag,
     4                pnch)
c
      implicit integer (a-z)
      real*8 basis, array
      real*8 y, w, wdiag
      real*8 hbar, omega, thresh, cnverg, eps
      real*8 energy, bigomg, lenscl, enscl
      real*8 sctlen, eta, scale, t0, gett0
      real*4 secnds, delta(20), time(20)
      character*(*) timh0, sptimv, i0stat, precon, card, coord
      character*80 title
      logical spac, prnt, hdiag, pnch
      dimension q(mgr,dim), qm(mgr,dim), wtm(mgr,dim)
      dimension pm(mgr,mgr,dim)    
      dimension nptt(mgr,dim), nmaxx(mgr,dim), ug(dim)
      dimension ham(mgr,dim), v(mgr,dim), trn(mgr,dim), eig(mgr,dim)
      dimension zq(4), zqm(4), zwtm(4), zpm(4)
      dimension zham(4), zv(4), ztrn(4), zeig(4)
      dimension nmax(4), npt(4), prnt(*)
      common/io/inp, iout      
      pointer(pbasis,basis(1))
      pointer(parray,array(1))
      pointer (py,y(1)), (py,iy(1))
      pointer (pw,w(1)), (pw,iw(1))
      pointer (pdiag,wdiag(1))
c
c
c        set pointers to the main grid used
c
      do 10 i=1,dim
         npt(i)=nptt(ug(i),i)
         nmax(i)=nmaxx(ug(i),i)
         zq(i) = q(ug(i),i)
         zqm(i) = qm(ug(i),i)
         zwtm(i) = wtm(ug(i),i)
         zpm(i) = pm(ug(i),ug(i),i)
         zham(i) = ham(ug(i),i)  
         zv(i) = v(ug(i),i)  
 10   continue
      n3d=1
      do 20 i=2,dim
         ztrn(i) = trn(ug(i),i)   
         zeig(i) = eig(ug(i),i)   
         n3d=n3d*nmax(i)
 20   continue
      ntot=n3d*nmax(1)   
      m=2*ntot
      if(itsolv) then
         maxit=intkey(card,'maximum-number-of-iterations',m,' ')
         maxvec=intkey(card,'maximum-number-of-vectors',m,' ')
      endif
      nrhs=1
      voft=1
      psi0=voft+ntot 
      driver=psi0+2*n3d
      need=wpadti(driver+m)
      call memory(need,py,ny,'pot',0)
c
c     form the one-body hamiltonian
c
      time(1)=secnds(0.0)
      if(tim.eq.1) then
         do 30 i=2,dim
            call addd(array(zham(i)),array(zham(i)),
     1                array(zv(i)),nmax(i))
 30      continue
      endif
      time(2)=secnds(0.0)
      delta(1)=time(2)-time(1)
c
c----------------------------------------------------------------------c
c     Calculate the time dependent perturbation.  It consists of a     c
c     space and a time part.                                           c
c----------------------------------------------------------------------c
      call vpert(y(voft),array(zv(1)),basis(zqm(2)),basis(zqm(3)),
     1           basis(zqm(4)),basis(zqm(1)),n3d,nmax(2),nmax(1),
     2           spdim,sptimv,timh0)
      time(3)=secnds(0.0)
      delta(2)=time(3)-time(2)
c----------------------------------------------------------------------c
c                                                                      
      t0=gett0(basis(zq(1)),npt(1))
c
c     calculate the spatial part of the incident wavefunction
c
      call mkpsi0(array(ztrn(2)),array(ztrn(3)),array(ztrn(4)),
     1            array(zeig(2)),array(zeig(3)),array(zeig(4)),
     2            basis(zpm(2)),basis(zpm(3)),basis(zpm(4)),
     3            basis(zqm(2)),basis(zqm(3)),basis(zqm(4)),
     4            basis(zwtm(2)),basis(zwtm(3)),basis(zwtm(4)),
     5            y(psi0),energy,t0,n3d,nmax(2),state,
     6            spdim,coord,tim,i0stat,prnt(4))
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
      call mkhpsi(array(zham(2)),array(zham(3)),array(zham(4)),
     1            y(driver),y(psi0),y(voft),basis(zpm(1)),
     2            basis(zwtm(1)),energy,n3d,nmax(2),nmax(3),
     3            nmax(4),nmax(1),spdim,i0stat,prnt(5))
      time(5)=secnds(0.0)
      delta(4)=time(5)-time(4)
      if(.not.itsolv) then
         ind=1
         hamil=iadtwp(ind+ntot)
         jpvt=wpadti(hamil+m*m)
         need=jpvt+m
         call memory(need,pw,nw,'direct-solve',0)
         call hamtot(w(hamil),array(zham(2)),array(zham(3)),
     1               array(zham(4)),array(zham(1)),y(voft),
     2               iw(ind),nmax(2),nmax(3),nmax(4),nmax(1),
     3               ntot,m,spdim,prnt(6))
         time(6)=secnds(0.0)
         delta(5)=time(6)-time(5)
         if(hdiag) then
            eigr=1
            eigi=eigr+m
            vecl=eigi+m
            vecr=vecl+m*m
            need=wpadti(vecr+m*m)
            call memory(need,pdiag,nwdiag,'diag',0)
            call eigslv(w(hamil),wdiag(eigr),wdiag(eigi),wdiag(vecl),
     1                  wdiag(vecr),m)
            time(7)=secnds(0.0)
            delta(6)=time(7)-time(6)
            call memory(-nwdiag,pdiag,idum,'diag',idum)
         else
            call cmpslv(w(hamil),y(driver),iw(jpvt),
     1                  m,1,prnt(7))
            time(7)=secnds(0.0)
            delta(6)=time(7)-time(6)
            call pltwfn(y(driver),basis(zqm(2)),basis(zqm(3)),
     1                  basis(zqm(4)),basis(zqm(1)),basis(zpm(2)),
     2                  basis(zpm(3)),basis(zpm(4)),basis(zpm(1)),
     3                  y(psi0),t0,energy,spdim,ntot,n3d,nmax(2),
     4                  nmax(3),nmax(4),nmax(1),timh0,i0stat,pnch)
            call initpsi(y(driver),y(psi0),basis(zpm(1)),n3d,nmax(1))
         endif
         call memory(-nw,pw,iidum,'direct-solve',iidum)
      else
         begin=1
         if(precon.ne.'none') then
            if(nblck.ge.ntot) then
               nblck=ntot
            endif
            hsub=1
            ind=wpadti(hsub+2*nblck*nblck)
            jnd=ind+(spdim+1)*ntot
            begin=iadtwp(jnd+ntot)
         endif
         pvec=begin
         vtrial=begin
         hpvec=pvec+m*maxvec
         tstor=hpvec
         h=hpvec+m*maxvec
         htmp=h+maxvec*maxvec
         b=htmp+maxvec*maxvec
         btmp=b+nrhs*maxvec
         resid=btmp+nrhs*maxvec
         soln=resid+m*maxvec
         t1=soln+m*nrhs
         t2=t1+m*maxvec
         need=wpadti(t2+m*maxvec)
         call memory(need,pw,nw,'iterative-solve',0)
         write(iout,*) precon
         if(precon.eq.'block'.or.precon.eq.'forced-block') then
            call cfctrs(w(hsub),y(driver),w(vtrial),w(tstor),
     1                  array(zham(2)),array(zham(3)),array(zham(4)),
     2                  array(zham(1)),y(voft),iw(ind),iw(jnd),
     3                  nmax(2),nmax(3),nmax(4),nmax(1),
     4                  ntot,nblck,spdim,1,prnt(8),.true.)
            time(6)=secnds(0.0)
            delta(5)=time(6)-time(5)
            ntrial=1 
         else
            call trials(w(pvec),array(zeig(2)),array(zeig(3)),
     1                  array(zeig(4)),array(ztrn(2)),array(ztrn(3)),
     2                  array(ztrn(4)),y(driver),thresh,nmax(2),nmax(3),
     3                  nmax(4),nmax(1),ntrial,m,maxvec,
     4                  spdim,.false.)
            time(6)=secnds(0.0)
            delta(5)=time(6)-time(5)
         endif
         scale=1.d0
         call lindvr(array(zham(2)),array(zham(3)),array(zham(4)),
     1               array(zham(1)),array(zeig(2)),array(zeig(3)),
     2               array(zeig(4)),array(ztrn(2)),array(ztrn(3)),
     3               array(ztrn(4)),y(voft),w(pvec),w(hpvec),w(h),
     4               w(htmp),w(b),w(btmp),y(driver),w(pvec),w(resid),
     5               w(soln),w(t1),w(t2),w(hsub),iw(ind),scale,
     6               cnverg,thresh,eps,nmax(2),nmax(3),nmax(4),
     7               nmax(1),spdim,m,ntot,n3d,nrhs,ntrial,maxit,
     8               maxvec,prnt(9),precon)
         time(7)=secnds(0.0)
         delta(6)=time(7)-time(6)
         call pltwfn(y(driver),basis(zqm(2)),basis(zqm(3)),
     1               basis(zqm(4)),basis(zqm(1)),basis(zpm(2)),
     2               basis(zpm(3)),basis(zpm(4)),basis(zpm(1)),
     3               y(psi0),t0,energy,spdim,ntot,n3d,nmax(2),
     4               nmax(3),nmax(4),nmax(1),timh0,i0stat,pnch)
         call initpsi(y(driver),y(psi0),basis(zpm(1)),n3d,nmax(1))
         call memory(-nw,pw,iidum,'iterative-solve',iidum)
      endif
c
      call memory(-ny,py,iidum,'pot',iidum)
      write(iout,1) (delta(i),i=1,4)
      if(.not.itsolv) then
         write(iout,2) delta(5), delta(6)
      else
         write(iout,3) delta(5), delta(6)
      endif
      return
 1    format(/,1x,'time for one body operator = ',e15.8,
     1       /,1x,'time for perturbation      = ',e15.8,
     2       /,1x,'time for right hand side   = ',e15.8,
     3       /,1x,'time for driving term      = ',e15.8)
 2    format(/,1x,'time for hamiltonian formation      = ',e15.8,
     1       /,1x,'time for direct eigen or linear system solve = ',
     2                                                     e15.8)
 3    format(/,1x,'time for trials or block LU factorization = ',e15.8,
     1       /,1x,'time iterative linear system solve        = ',e15.8)
      end
