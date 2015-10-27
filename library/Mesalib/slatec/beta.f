*deck beta
      function beta (a, b)
c***begin prologue  beta
c***purpose  compute the complete beta function.
c***library   slatec (fnlib)
c***category  c7b
c***type      single precision (beta-s, dbeta-d, cbeta-c)
c***keywords  complete beta function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c beta computes the complete beta function.
c
c input parameters:
c       a   real and positive
c       b   real and positive
c
c***references  (none)
c***routines called  albeta, gamlim, gamma, r1mach, xermsg
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900727  added external statement.  (wrb)
c***end prologue  beta
      external gamma
      save xmax, alnsml
      data xmax, alnsml /0., 0./
c***first executable statement  beta
      if (alnsml.ne.0.0) go to 10
      call gamlim (xmin, xmax)
      alnsml = log(r1mach(1))
c
 10   if (a .le. 0. .or. b .le. 0.) call xermsg ('slatec', 'beta',
     +   'both arguments must be gt 0', 2, 2)
c
      if (a+b.lt.xmax) beta = gamma(a) * gamma(b) / gamma(a+b)
      if (a+b.lt.xmax) return
c
      beta = albeta (a, b)
      if (beta .lt. alnsml) call xermsg ('slatec', 'beta',
     +   'a and/or b so big beta underflows', 1, 2)
c
      beta = exp (beta)
c
      return
      end
