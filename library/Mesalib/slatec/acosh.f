*deck acosh
      function acosh (x)
c***begin prologue  acosh
c***purpose  compute the arc hyperbolic cosine.
c***library   slatec (fnlib)
c***category  c4c
c***type      single precision (acosh-s, dacosh-d, cacosh-c)
c***keywords  acosh, arc hyperbolic cosine, elementary functions, fnlib,
c             inverse hyperbolic cosine
c***author  fullerton, w., (lanl)
c***description
c
c acosh(x) computes the arc hyperbolic cosine of x.
c
c***references  (none)
c***routines called  r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  acosh
      save aln2,xmax
      data aln2 / 0.6931471805 5994530942e0/
      data xmax /0./
c***first executable statement  acosh
      if (xmax.eq.0.) xmax = 1.0/sqrt(r1mach(3))
c
      if (x .lt. 1.0) call xermsg ('slatec', 'acosh', 'x less than 1',
     +   1, 2)
c
      if (x.lt.xmax) acosh = log (x + sqrt(x*x-1.0))
      if (x.ge.xmax) acosh = aln2 + log(x)
c
      return
      end
