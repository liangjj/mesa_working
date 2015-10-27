*deck sdatrp
      subroutine sdatrp (x, xout, yout, ypout, neq, kold, phi, psi)
c***begin prologue  sdatrp
c***subsidiary
c***purpose  interpolation routine for sdassl.
c***library   slatec (dassl)
c***type      single precision (sdatrp-s, ddatrp-d)
c***author  petzold, linda r., (llnl)
c***description
c-----------------------------------------------------------------------
c     the methods in subroutine sdastp use polynomials
c     to approximate the solution. sdatrp approximates the
c     solution and its derivative at time xout by evaluating
c     one of these polynomials, and its derivative,there.
c     information defining this polynomial is passed from
c     sdastp, so sdatrp cannot be used alone.
c
c     the parameters are:
c     x     the current time in the integration.
c     xout  the time at which the solution is desired
c     yout  the interpolated approximation to y at xout
c           (this is output)
c     ypout the interpolated approximation to yprime at xout
c           (this is output)
c     neq   number of equations
c     kold  order used on last successful step
c     phi   array of scaled divided differences of y
c     psi   array of past stepsize history
c-----------------------------------------------------------------------
c***routines called  (none)
c***revision history  (yymmdd)
c   830315  date written
c   901009  finished conversion to slatec 4.0 format (f.n.fritsch)
c   901019  merged changes made by c. ulrich with slatec 4.0 format.
c   901026  added explicit declarations for all variables and minor
c           cosmetic changes to prologue.  (fnf)
c***end prologue  sdatrp
c
      integer  neq, kold
      real  x, xout, yout(*), ypout(*), phi(neq,*), psi(*)
c
      integer  i, j, koldp1
      real  c, d, gamma, temp1
c
c***first executable statement  sdatrp
      koldp1=kold+1
      temp1=xout-x
      do 10 i=1,neq
         yout(i)=phi(i,1)
10       ypout(i)=0.0e0
      c=1.0e0
      d=0.0e0
      gamma=temp1/psi(1)
      do 30 j=2,koldp1
         d=d*gamma+c/psi(j-1)
         c=c*gamma
         gamma=(temp1+psi(j-1))/psi(j)
         do 20 i=1,neq
            yout(i)=yout(i)+c*phi(i,j)
20          ypout(i)=ypout(i)+d*phi(i,j)
30       continue
      return
c
c------end of subroutine sdatrp------
      end
