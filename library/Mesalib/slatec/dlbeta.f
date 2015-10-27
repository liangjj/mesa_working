*deck dlbeta
      double precision function dlbeta (a, b)
c***begin prologue  dlbeta
c***purpose  compute the natural logarithm of the complete beta
c            function.
c***library   slatec (fnlib)
c***category  c7b
c***type      double precision (albeta-s, dlbeta-d, clbeta-c)
c***keywords  fnlib, logarithm of the complete beta function,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c dlbeta(a,b) calculates the double precision natural logarithm of
c the complete beta function for double precision arguments
c a and b.
c
c***references  (none)
c***routines called  d9lgmc, dgamma, dlngam, dlnrel, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900727  added external statement.  (wrb)
c***end prologue  dlbeta
      double precision a, b, p, q, corr, sq2pil, d9lgmc, dgamma, dlngam,
     1  dlnrel
      external dgamma
      save sq2pil
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
c***first executable statement  dlbeta
      p = min (a, b)
      q = max (a, b)
c
      if (p .le. 0.d0) call xermsg ('slatec', 'dlbeta',
     +   'both arguments must be gt zero', 1, 2)
c
      if (p.ge.10.d0) go to 30
      if (q.ge.10.d0) go to 20
c
c p and q are small.
c
      dlbeta = log (dgamma(p) * (dgamma(q)/dgamma(p+q)) )
      return
c
c p is small, but q is big.
c
 20   corr = d9lgmc(q) - d9lgmc(p+q)
      dlbeta = dlngam(p) + corr + p - p*log(p+q)
     1  + (q-0.5d0)*dlnrel(-p/(p+q))
      return
c
c p and q are big.
c
 30   corr = d9lgmc(p) + d9lgmc(q) - d9lgmc(p+q)
      dlbeta = -0.5d0*log(q) + sq2pil + corr + (p-0.5d0)*log(p/(p+q))
     1  + q*dlnrel(-p/(p+q))
      return
c
      end
