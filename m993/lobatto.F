*deck lobatto.f 
c***begin prologue     lobatto
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           lobatto functions
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            1. calculate piecewise lobatto dvr functions and
c***                      their one-body matrices
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       lobatto
      subroutine lobatto(phamil,mass,edge,coord,vkey,qtyp,typpot,
     1                   bcleft,bcright,n,npt,nrq,nreg,
     2                   nodiag,nword,prn,pgrid,ngot)
c
      implicit integer (a-z)
      character*(*) coord, vkey, qtyp, typpot
      character*24 typarg
      real*8 grid, mat, glbg, ham, hamil
      real*8 edge, mass, dscale  
#ifdef DECPOINTER
      integer*8 phamil
      integer*8 pgrid
#endif DECPOINTER
#ifdef SGIPOINTER
      integer*4 phamil
      integer*4 pgrid
#endif SGIPOINTER
      logical prn, nodiag
      dimension pgrid(nreg,2), npt(nreg), nrq(nreg)
      dimension ngot(nreg,2), edge(nreg+1), prn(*)
      dimension words(2)
      common/io/inp, iout      
      pointer (phamil,hamil(1))
      pointer (pglbg,glbg(1))
      pointer (pham,ham(1))
      pointer (pg,grid(1))
      pointer (pm,mat(1))
      typarg='linear'
      if(coord.eq.'cylindrical') then
         typarg='quadratic'
      endif 
      n=0
      do 10 i=1,nreg
c
c        number of internal functions
c
         n = n + npt(i) - 2
 10   continue   
c
c     add one bridge function between each interval.
c
      n = n + nreg - 1 
c
c     add the extreme left and extreme right points
c     we have not yet dropped any functions at the endpoints.
c
      n = n + 2
c
c     get memory for the global arrays
c
      x=1
      xwt=x+n
      px=xwt+n
      dpx=px+n*n
      ddpx=dpx+n*n                 
      norm=ddpx+n*n
      wrd=norm+n-1
      need=wpadti(wrd+1)
      call getmem(need,pglbg,words(1),'global',0)
      call rzero(glbg,wrd)
      do 20 i=1,nreg
         q=1
         wt=q+npt(i)
         p=wt+npt(i)
         dp=p+npt(i)*npt(i)
         ddp=dp+npt(i)*npt(i)
         need=wpadti(ddp+npt(i)*npt(i)) 
         call getmem(need,pgrid(i,1),ngot(i,1),'grid',0)
         pg=pgrid(i,1)
c
c        calculate the un-normalized sector functions and their
c        derivatives.
c
         call drvply(grid(q),grid(wt),grid(p),grid(dp),
     1               grid(ddp),edge(i),qtyp,coord,npt(i)-1,
     2               nrq(i),prn(1))

c        calculate the overlap, bloch and kinetic energy matrix
c        over un-normalized functions.

         ov=1
         bl=ov+npt(i)*npt(i)
         ke=bl+npt(i)*npt(i)                 
         need=wpadti(ke+npt(i)*npt(i))                 
         call getmem(need,pgrid(i,2),ngot(i,2),'mat',0)
         pm=pgrid(i,2)
         call ovmat(mat(ov),grid(p),grid(p),grid(wt),
     #                   npt(i),prn(3),i)
         call blmat(mat(bl),grid(p),grid(dp),grid(q),
     #              mass,typarg,npt(i),prn(3),i)
         call kemat(mat(ke),mat(bl),grid(p),grid(dp),grid(ddp),
     #              grid(q),grid(wt),mass,typarg,npt(i),prn(3),i)
c
c
c
 20   continue   
c
c     calculate the global functions
c
c----------------------------------------------------------------------c
c
c                  x = points  xwt=weights px=polynomials
c                  dpx=first derivatives  ddpx=second derivatives
c
c                  these are all full.  no removal of endpoints.
c----------------------------------------------------------------------c
      call conbsis(pgrid(1,1),glbg(x),glbg(xwt),glbg(px),glbg(dpx),
     1             glbg(ddpx),glbg(norm),edge,n,npt,nreg,prn(4))
c
c     calculate the matrix elements and potentials
c
      hmat=1
      v=hmat+n*n
      need=wpadti(v+n)
      call getmem(need,pham,words(2),'ham',0)
      call vone(vkey,glbg(x),ham(v),dscale,typpot,n,.true.,prn(6))
c----------------------------------------------------------------------c
c                  hmat=hamiltonian v=potential
c                  first and last functions not removed.
c----------------------------------------------------------------------c
      call conham(pgrid(1,2),ham(hmat),ham(v),glbg(norm),n,npt,
     #            nreg,prn(7))
c
c     now make all of the physical quantites and diagonalize the
c                          hamiltonian.
c
      ntot=n
      start=1
      end=n
      if(bcleft.eq.0) then
          ntot=ntot-1
          start=start+1
      endif
      if(bcright.eq.0) then
         ntot=ntot-1
         end=end-1
      endif 
      h=1
      vphy=h+ntot*ntot
      h0=vphy+ntot
      srf=h0+ntot*ntot
      q0=srf+2
      q1=q0+1
      qwt=q1+ntot
      pq=qwt+ntot
      dpq=pq+ntot*ntot
      ddpq=dpq+ntot*ntot            
      need=ddpq+ntot*ntot
      if(.not.nodiag) then
         eigv=need
         rgama=eigv+ntot*ntot
         eig=rgama+ntot*2
         eigv0=eig+ntot
         rgama0=eigv0+ntot*ntot
         eig0=rgama0+2*ntot     
         junk=eig0+ntot
         need=junk+5*ntot
      endif
      need=wpadti(need)
      call getmem(need,phamil,nword,'diag',0) 
      call hamphy(glbg(x),glbg(xwt),glbg(px),glbg(dpx),glbg(ddpx),
     1            ham(hmat),ham(v),hamil(q0),hamil(q1),hamil(qwt),
     2            hamil(pq),hamil(dpq),hamil(ddpq),hamil(h),hamil(vphy),
     3            hamil(srf),n,ntot,start)
c
      if(.not.nodiag) then
         call diarep(hamil(h),hamil(eig),hamil(eigv),
     1               hamil(junk),ntot,prn(9))
         call rmamp(hamil(eigv),hamil(srf),hamil(rgama),ntot,prn(11))
      endif
      call copy(hamil(h),hamil(h0),ntot*ntot)
      call vec2di(hamil(h0),hamil(h0),hamil(vphy),'subtract',ntot)
      if(.not.nodiag) then
         call diarep(hamil(h0),hamil(eig0),hamil(eigv0),hamil(junk),
     1               ntot,prn(9))
         call rmamp(hamil(eigv0),hamil(srf),hamil(rgama0),ntot,prn(11))
      endif
c
c                
      do 30 i=1,nreg
         call getmem(-ngot(i,1),pgrid(i,1),idum,'grid',idum)
         call getmem(-ngot(i,2),pgrid(i,2),idum,'mat',idum)
 30   continue   
      call getmem(-words(1),pglbg,dum,'global',dum)
      call getmem(-words(2),pham,dum,'ham',dum)
      return
      end
















