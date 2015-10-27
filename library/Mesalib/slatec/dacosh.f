*deck dacosh
      double precision function dacosh (x)
c***begin prologue  dacosh
c***purpose  compute the arc hyperbolic cosine.
c***library   slatec (fnlib)
c***category  c4c
c***type      double precision (acosh-s, dacosh-d, cacosh-c)
c***keywords  acosh, arc hyperbolic cosine, elementary functions, fnlib,
c             inverse hyperbolic cosine
c***author  fullerton, w., (lanl)
c***description
c
c dacosh(x) calculates the double precision arc hyperbolic cosine for
c double precision argument x.  the result is returned on the
c positive branch.
c
c***references  (none)
c***routines called  d1mach, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dacosh
      double precision x, dln2, xmax,  d1mach
      save dln2, xmax
      data dln2 / 0.6931471805 5994530941 7232121458 18 d0 /
      data xmax / 0.d0 /
c***first executable statement  dacosh
      if (xmax.eq.0.d0) xmax = 1.0d0/sqrt(d1mach(3))
c
      if (x .lt. 1.d0) call xermsg ('slatec', 'dacosh',
     +   'x less than 1', 1, 2)
c
      if (x.lt.xmax) dacosh = log (x+sqrt(x*x-1.0d0))
      if (x.ge.xmax) dacosh = dln2 + log(x)
c
      return
      end
