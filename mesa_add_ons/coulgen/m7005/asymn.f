*deck asymn
      subroutine asymn(gl,dgl,e0,r,rinv,eta,l,npt,n,level)
c***begin prologue     asymn
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            asymntotic expansion for irregular negative
c***                   energy coulomb function at large rho.
c***
c***description        exponentially decaying coulomb function is
c***                   computed using an asymptotic expansion at
c***                   large rho.
c***references
c
c***routines called
c
c***end prologue       asymn
c
      implicit integer (a-z)
      real*8 gl, dgl, e0, eta, tol, r, rinv, pre, add, fac
      real*8 zero, half, one, dsqrt2
      parameter (tol=1.d-14)
      dimension gl(npt), dgl(npt), e0(0:n), r(npt), rinv(npt)
      common/io/inp,iout
      data half, one, dsqrt2 / .5d0, 1.d0, .70710678118654752440d0/
      call rzero(gl,npt)
      call rzero(dgl,npt)
      do 10 pt=1,npt
c**********************************************************************c
c            gl = -1./sqrt(2.) * (series) * exp(-rho)                  c
c                              * (2.rho)**(-eta)                       c
c                                                                      c
c            series = sum ( n =0 to n= infinity ) e0(n)                c
c                             * (-2.*rho)**(-eta)                      c
c                                                                      c
c                * Henry incorrectly uses 2.*rho                       c
c**********************************************************************c
c**********************************************************************c
c                do the series expansion to the level tol              c
c**********************************************************************c
         pre=one
         fac=-half*rinv(pt)
         do 20 i=0,n
            add=e0(i)*pre
            if (abs(add).le.tol) go to 30
            gl(pt)=gl(pt)+add
            dgl(pt)=dgl(pt)+i*add
            pre=pre*fac
   20    continue
   30    if (level.ge.3) then
             write(iout,*) 'asymptotic expansion converged in ',i,
     1                     ' terms'
         endif
   10 continue     
      do 40 pt=1,npt
         fac=half*rinv(pt)
         pre=dsqrt2*exp(-r(pt))
         add=fac**eta
         gl(pt)=-gl(pt)*pre*add
         dgl(pt)=-gl(pt)*(one+eta*rinv(pt))-pre*add*dgl(pt)*rinv(pt)
   40 continue
      return
      end

