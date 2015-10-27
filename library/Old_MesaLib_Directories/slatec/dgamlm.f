*deck dgamlm
      subroutine dgamlm (xmin, xmax)
c***begin prologue  dgamlm
c***purpose  compute the minimum and maximum bounds for the argument in
c            the gamma function.
c***library   slatec (fnlib)
c***category  c7a, r2
c***type      double precision (gamlim-s, dgamlm-d)
c***keywords  complete gamma function, fnlib, limits, special functions
c***author  fullerton, w., (lanl)
c***description
c
c calculate the minimum and maximum legal bounds for x in gamma(x).
c xmin and xmax are not the only bounds, but they are the only non-
c trivial ones to calculate.
c
c             output arguments --
c xmin   double precision minimum legal value of x in gamma(x).  any
c        smaller value of x might result in underflow.
c xmax   double precision maximum legal value of x in gamma(x).  any
c        larger value of x might cause overflow.
c
c***references  (none)
c***routines called  d1mach, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dgamlm
      double precision xmin, xmax, alnbig, alnsml, xln, xold, d1mach
c***first executable statement  dgamlm
      alnsml = log(d1mach(1))
      xmin = -alnsml
      do 10 i=1,10
        xold = xmin
        xln = log(xmin)
        xmin = xmin - xmin*((xmin+0.5d0)*xln - xmin - 0.2258d0 + alnsml)
     1    / (xmin*xln+0.5d0)
        if (abs(xmin-xold).lt.0.005d0) go to 20
 10   continue
      call xermsg ('slatec', 'dgamlm', 'unable to find xmin', 1, 2)
c
 20   xmin = -xmin + 0.01d0
c
      alnbig = log (d1mach(2))
      xmax = alnbig
      do 30 i=1,10
        xold = xmax
        xln = log(xmax)
        xmax = xmax - xmax*((xmax-0.5d0)*xln - xmax + 0.9189d0 - alnbig)
     1    / (xmax*xln-0.5d0)
        if (abs(xmax-xold).lt.0.005d0) go to 40
 30   continue
      call xermsg ('slatec', 'dgamlm', 'unable to find xmax', 2, 2)
c
 40   xmax = xmax - 0.01d0
      xmin = max (xmin, -xmax+1.d0)
c
      return
      end
