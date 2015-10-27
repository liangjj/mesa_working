      subroutine df(t,psi,dpsi,dum,idum,neq,h01,h02,h03,eig1,eig2,eig3,
     1              v0,vt,tpert,efield,omega,width,shift,hbar,
     2              dim,nmax,n) 
      implicit integer(a-z)
      real*8 t, psi, dpsi, dum, ddum
      real*8 h01, h02, h03, eig1, eig2, eig3
      real*8 v0, vt, omega, efield, width, shift
      real*8 hbar
      character*(*) tpert
      dimension nmax(3)
      dimension h01(nmax(1),*), h02(nmax(2),*), h03(nmax(3),*)
      dimension eig1(nmax(1)), eig2(nmax(2)), eig3(nmax(3))
      dimension v0(n), vt(n)
      dimension psi(neq), dpsi(neq)
      common/io/inp, iout
c
c     psi( dpsi ) is to be regarded as a vector containing the real
c     ( imaginary ) part of psi ( dpsi ) in the first ( last ) n locations.
c     it can be thought of as an array of the form array(n,2)
c     where
      n=neq/2
c
c     get the time dependent part of the potential and add it to the
c     time independent part.
      call vpert(vt,eig1,eig2,eig3,t,efield,omega,width,shift,
     1           nmax,dim,tpert)
      call vadd(vt,vt,v0,n)
c
c     calculate derivative
c
      if(dim.eq.1) then
         call honv(h01,ddum,ddum,vt,psi(n+1),dpsi,n,nmax(1),
     1             iidum,iidum,dim)
         call honv(h01,ddum,ddum,vt,psi,dpsi(n+1),n,nmax(1),
     1             iidum,iidum,dim)
      elseif(dim.eq.2) then
         call honv(h02,h01,ddum,vt,psi(n+1),dpsi,n,nmax(2),
     1             nmax(1),iidum,dim)
         call honv(h02,h01,ddum,vt,psi,dpsi(n+1),n,nmax(2),
     1             nmax(1),iidum,dim)
      elseif(dim.eq.3) then
         call honv(h03,h02,h01,vt,psi(n+1),dpsi,n,nmax(3),
     1             nmax(2),nmax(1),dim)
         call honv(h03,h02,h01,vt,psi,dpsi(n+1),n,nmax(3),
     1             nmax(2),nmax(1),dim)
      endif
      call vscale(dpsi,dpsi,1.d0/hbar,n)
      call vscale(dpsi(n+1),dpsi(n+1),-1.d0/hbar,n)
      return
      end
