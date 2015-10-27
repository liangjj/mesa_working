*deck psit.f 
c***begin prologue     psit
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time, dvr, orthogonal polynomial
c***                   bose-einstein
c***author             schneider, b. i.(nsf)
c***source             tprop
c***purpose            linear time dependent schroedinger equation 
c***                   using time based dvr.
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       psit
      subroutine psit(pham,pvt,spdim,dim,nc,nphy,key,typpot,spac,tim,
     1                i0stat,state,coord,itsolv,cnverg,thresh,
     2                eps,precon,nblck,nwpot,card,prnt,hdiag,pnch)
c
      implicit integer (a-z)
      integer*8 pham, pvt, pv
      real*8 ham1, ham2, ham3, ham4, v, vt, y, w, wdiag
      real*8 thresh, cnverg, eps
      real*8 energy
      real*8 scale, t0
      real*4 secnds, delta(20), time(20)
      character*(*) i0stat, key, typpot, precon, card, coord
      character*80 title
      logical spac, prnt, hdiag, pnch
      dimension pham(dim)
      dimension h0(4), eigv0(4), eig0(4)
      dimension q0(4), q(4), qwt(4)
      dimension pq(4), dpq(4), nphy(4), prnt(*), key(4), typpot(4)
      common/io/inp, iout     
      pointer(ph1,ham1(1))
      pointer(ph2,ham2(1))
      pointer(ph3,ham3(1))
      pointer(ph4,ham4(1))
      pointer (pv,v(1))
      pointer (pvt,vt(1))
      pointer (py,y(1)), (py,iy(1))
      pointer (pw,w(1)), (pw,iw(1))
      pointer (pdiag,wdiag(1))
c
c     get the pure time dependent part of the hamiltonian.
c
      call hamt(pham(dim),pvt,key(dim),typpot(dim),
     1          nphy(dim),nwpot,prn)
c
c     set the local pointers to their global value
c
      do 10 i=1,dim
         notim=.true.
         if(i.eq.dim) then
            notim=.false.
         endif
         call set(pham(i),pvt,h0(i),eigv0(i),eig0(i),q0(i),
     1            q(i),qwt(i),pq(i),dpq(i),nphy(i),prnt(1),notim,spac)
 10   continue   
      ph1=pham(1)
      if(spdim.gt.1) then
         ph2=pham(2)
      endif
      if(spdim.gt.2) then
         ph3=pham(3)
      endif
      ph4=pham(dim)
c
      n3d=1
      do 20 i=1,spdim
         n3d=n3d*nphy(i)
 20   continue
      ntot=n3d*nphy(dim)   
      m=2*ntot
      if(itsolv) then
         maxit=intkey(card,'maximum-number-of-iterations',m,' ')
         maxvec=intkey(card,'maximum-number-of-vectors',m,' ')
      endif
      nrhs=1
      psi0=1 
      driver=psi0+2*n3d*nc
      need=wpadti(driver+m*nc)
      call memory(need,py,ny,'pot',0)
c                                                                      
      t0=ham4(q0(dim))
c
c     calculate the spatial part of the incident wavefunction
c
      time(1)=secnds(0.0)
      call mkpsi0(ham1(eigv0(1)),ham2(eigv0(2)),ham3(eigv0(3)),
     1            ham1(eig0(1)),ham2(eig0(2)),ham3(eig0(3)),
     2            ham1(pq(1)),ham2(pq(2)),ham3(pq(3)),
     3            ham1(q(1)),ham2(q(2)),ham3(q(3)),
     4            ham1(qwt(1)),ham2(qwt(2)),ham3(qwt(3)),
     5            y(psi0),energy,t0,n3d,nphy(1),nc,
     6            state,spdim,coord,tim,i0stat,prnt(4))
      time(2)=secnds(0.0)
      delta(1)=time(2)-time(1)
c
c----------------------------------------------------------------------c
c     Calculate the time dependent perturbation.  It consists of a     c
c     space and a time part.                                           c
c----------------------------------------------------------------------c
      call vpert(pv,pvt,ham1(q(1)),ham2(q(2)),ham3(q(3)),
     1           ham4(q(dim)),ham1(pq(1)),ham2(pq(2)),ham3(pq(3)),
     2           y(driver),n3d,nphy(1),nphy(dim),nc,spdim,coord,tim)
      time(3)=secnds(0.0)
      delta(2)=time(3)-time(2)
c----------------------------------------------------------------------c
      call mkhpsi(ham1(h0(1)),ham2(h0(2)),ham3(h0(3)),
     1            y(driver),y(psi0),v,ham4(pq(dim)),
     2            ham4(qwt(dim)),energy,n3d,nphy(1),nphy(2),nphy(3),
     3            nphy(dim),spdim,prnt(5))
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
      if(.not.itsolv) then
         ind=1
         hamil=iadtwp(ind+ntot)
         jpvt=wpadti(hamil+m*m)
         need=jpvt+m
         call memory(need,pw,nw,'direct-solve',0)
         call hamtot(w(hamil),ham1(h0(1)),ham2(h0(2)),
     1               ham3(h0(3)),ham4(h0(dim)),v,
     2               iw(ind),nphy(1),nphy(2),nphy(3),nphy(dim),
     3               ntot,m,spdim,prnt(6))
         time(5)=secnds(0.0)
         delta(4)=time(5)-time(4)
         if(hdiag) then
            eigr=1
            eigi=eigr+m
            vecl=eigi+m
            vecr=vecl+m*m
            need=wpadti(vecr+m*m)
            call memory(need,pdiag,nwdiag,'diag',0)
            call eigslv(w(hamil),wdiag(eigr),wdiag(eigi),wdiag(vecl),
     1                  wdiag(vecr),m)
            time(6)=secnds(0.0)
            delta(5)=time(6)-time(5)
            call memory(-nwdiag,pdiag,idum,'diag',idum)
         else
            call cmpslv(w(hamil),y(driver),iw(jpvt),
     1                  m,1,prnt(7))
            time(6)=secnds(0.0)
            delta(5)=time(6)-time(5)
            call pltwfn(y(driver),ham1(q(1)),ham2(q(2)),
     1                  ham3(q(3)),ham4(q(dim)),ham1(pq(1)),
     2                  ham2(pq(2)),ham3(pq(3)),ham4(pq(dim)),
     3                  y(psi0),t0,energy,spdim,ntot,n3d,nphy(1),
     4                  nphy(2),nphy(3),nphy(dim),typpot(dim),
     5                  i0stat,pnch)
            call initpsi(y(driver),y(psi0),ham4(pq(dim)),n3d,nphy(dim))
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
     1                  ham1(h0(1)),ham2(h0(2)),ham3(h0(3)),
     2                  ham4(h0(dim)),v,iw(ind),iw(jnd),
     3                  nphy(1),nphy(2),nphy(3),nphy(dim),
     4                  ntot,nblck,spdim,1,prnt(8),.true.)
            time(5)=secnds(0.0)
            delta(4)=time(5)-time(4)
            ntrial=1 
         else
            call trials(w(pvec),ham1(eig0(1)),ham2(eig0(2)),
     1                  ham3(eig0(3)),ham1(eigv0(1)),ham2(eigv0(2)),
     2                  ham3(eigv0(3)),y(driver),thresh,nphy(1),nphy(2),
     3                  nphy(3),nphy(dim),ntrial,m,maxvec,
     4                  spdim,.false.)
            time(5)=secnds(0.0)
            delta(4)=time(5)-time(4)
         endif
         scale=1.d0
         call lindvr(ham1(h0(1)),ham2(h0(2)),ham3(h0(3)),
     1               ham4(h0(dim)),ham1(eig0(1)),ham2(eig0(2)),
     2               ham3(eig0(3)),ham1(eigv0(1)),ham2(eigv0(2)),
     3               ham3(eigv0(3)),v,w(pvec),w(hpvec),w(h),
     4               w(htmp),w(b),w(btmp),y(driver),w(pvec),w(resid),
     5               w(soln),w(t1),w(t2),w(hsub),iw(ind),scale,
     6               cnverg,thresh,eps,nphy(1),nphy(2),nphy(3),
     7               nphy(dim),spdim,m,ntot,n3d,nrhs,ntrial,maxit,
     8               maxvec,prnt(9),precon)
         time(6)=secnds(0.0)
         delta(5)=time(6)-time(5)
         call pltwfn(y(driver),ham1(q(1)),ham2(q(2)),ham3(q(3)),
     1               ham4(q(dim)),ham1(pq(1)),ham2(pq(2)),
     2               ham3(pq(3)),ham4(pq(dim)),y(psi0),t0,
     3               energy,spdim,ntot,n3d,nphy(1),nphy(2),nphy(3),
     4               nphy(dim),typpot(dim),i0stat,pnch)
         call initpsi(y(driver),y(psi0),ham4(pq(dim)),n3d,nphy(dim))
         call memory(-nw,pw,iidum,'iterative-solve',iidum)
      endif
c
      call memory(-ny,py,iidum,'pot',iidum)
      write(iout,1) (delta(i),i=1,3)
      if(.not.itsolv) then
         write(iout,2) delta(4), delta(5)
      else
         write(iout,3) delta(4), delta(5)
      endif
      return
 1    format(/,1x,'time for right hand side   = ',e15.8,
     1       /,1x,'time for perturbation      = ',e15.8,
     2       /,1x,'time for driving term      = ',e15.8)
 2    format(/,1x,'time for hamiltonian formation      = ',e15.8,
     1       /,1x,'time for direct eigen or linear system solve = ',
     2                                                     e15.8)
 3    format(/,1x,'time for trials or block LU factorization = ',e15.8,
     1       /,1x,'time iterative linear system solve        = ',e15.8)
      end






















