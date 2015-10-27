*deck gamlim
      subroutine gamlim (xmin, xmax)
c***begin prologue  gamlim
c***purpose  compute the minimum and maximum bounds for the argument in
c            the gamma function.
c***library   slatec (fnlib)
c***category  c7a, r2
c***type      single precision (gamlim-s, dgamlm-d)
c***keywords  complete gamma function, fnlib, limits, special functions
c***author  fullerton, w., (lanl)
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
c
c***references  (none)
c***routines called  r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  gamlim
c***first executable statement  gamlim
      alnsml = log(r1mach(1))
      xmin = -alnsml
      do 10 i=1,10
        xold = xmin
        xln = log(xmin)
        xmin = xmin - xmin*((xmin+0.5)*xln - xmin - 0.2258 + alnsml)
     1    / (xmin*xln + 0.5)
        if (abs(xmin-xold).lt.0.005) go to 20
 10   continue
      call xermsg ('slatec', 'gamlim', 'unable to find xmin', 1, 2)
c
 20   xmin = -xmin + 0.01
c
      alnbig = log(r1mach(2))
      xmax = alnbig
      do 30 i=1,10
        xold = xmax
        xln = log(xmax)
        xmax = xmax - xmax*((xmax-0.5)*xln - xmax + 0.9189 - alnbig)
     1    / (xmax*xln - 0.5)
        if (abs(xmax-xold).lt.0.005) go to 40
 30   continue
      call xermsg ('slatec', 'gamlim', 'unable to find xmax', 2, 2)
c
 40   xmax = xmax - 0.01
      xmin = max (xmin, -xmax+1.)
c
      return
      end
