*deck albeta
      function albeta (a, b)
c***begin prologue  albeta
c***purpose  compute the natural logarithm of the complete beta
c            function.
c***library   slatec (fnlib)
c***category  c7b
c***type      single precision (albeta-s, dlbeta-d, clbeta-c)
c***keywords  fnlib, logarithm of the complete beta function,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c albeta computes the natural log of the complete beta function.
c
c input parameters:
c       a   real and positive
c       b   real and positive
c
c***references  (none)
c***routines called  alngam, alnrel, gamma, r9lgmc, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900727  added external statement.  (wrb)
c***end prologue  albeta
      external gamma
      save sq2pil
      data sq2pil / 0.9189385332 0467274 e0 /
c***first executable statement  albeta
      p = min (a, b)
      q = max (a, b)
c
      if (p .le. 0.0) call xermsg ('slatec', 'albeta',
     +   'both arguments must be gt zero', 1, 2)
      if (p.ge.10.0) go to 30
      if (q.ge.10.0) go to 20
c
c p and q are small.
c
      albeta = log(gamma(p) * (gamma(q)/gamma(p+q)) )
      return
c
c p is small, but q is big.
c
 20   corr = r9lgmc(q) - r9lgmc(p+q)
      albeta = alngam(p) + corr + p - p*log(p+q) +
     1  (q-0.5)*alnrel(-p/(p+q))
      return
c
c p and q are big.
c
 30   corr = r9lgmc(p) + r9lgmc(q) - r9lgmc(p+q)
      albeta = -0.5*log(q) + sq2pil + corr + (p-0.5)*log(p/(p+q))
     1  + q*alnrel(-p/(p+q))
      return
c
      end
