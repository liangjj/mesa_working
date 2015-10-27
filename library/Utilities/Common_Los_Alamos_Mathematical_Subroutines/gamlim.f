      subroutine gamlim(xmin,xmax)
c***begin prologue  gamlim
c***date written   770401   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c7a,r2
c***keywords  gamma function,limits,special function
c***author  fullerton, w., (lanl)
c***purpose  computes the minimum and maximum bounds for x in gamma(x).
c***description
c
c calculate the minimum and maximum legal bounds for x in gamma(x).
c xmin and xmax are not the only bounds, but they are the only non-
c trivial ones to calculate.
c
c             output arguments --
c xmin   minimum legal value of x in gamma(x).  any smaller value of
c        x might result in underflow.
c xmax   maximum legal value of x in gamma(x).  any larger value will
c        cause overflow.
c***references  (none)
c***routines called  r1mach,xerror
c***end prologue  gamlim
c***first executable statement  gamlim
      implicit real*8(a-h,o-z)
      alnsml = log(r1mach(1))
      xmin = -alnsml
      do 10 i=1,10
        xold = xmin
        xln = log(xmin)
        xmin = xmin - xmin*((xmin+0.5d0)*xln - xmin - 0.2258d0 + alnsml)
     1    / (xmin*xln + 0.5d0)
        if (abs(xmin-xold).lt.0.005d0) go to 20
 10   continue
      call lnkerr ( 'gamlim  unable to find xmin')
c
 20   xmin = -xmin + 0.01d0
c
      alnbig = log(r1mach(2))
      xmax = alnbig
      do 30 i=1,10
        xold = xmax
        xln = log(xmax)
        xmax = xmax - xmax*((xmax-0.5d0)*xln - xmax + 0.9189d0 - alnbig)
     1    / (xmax*xln - 0.5d0)
        if (abs(xmax-xold).lt.0.005d0) go to 40
 30   continue
      call lnkerr ( 'gamlim  unable to find xmax')
c
 40   xmax = xmax - 0.01d0
      xmin = max (xmin, -xmax+1.d0)
c
      return
      end
