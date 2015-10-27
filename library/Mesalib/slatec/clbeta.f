*deck clbeta
      complex function clbeta (a, b)
c***begin prologue  clbeta
c***purpose  compute the natural logarithm of the complete beta
c            function.
c***library   slatec (fnlib)
c***category  c7b
c***type      complex (albeta-s, dlbeta-d, clbeta-c)
c***keywords  fnlib, logarithm of the complete beta function,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c clbeta computes the natural log of the complex valued complete beta
c function of complex parameters a and b.  this is a preliminary version
c which is not accurate.
c
c input parameters:
c       a   complex and the real part of a positive
c       b   complex and the real part of b positive
c
c***references  (none)
c***routines called  clngam, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  clbeta
      complex a, b, clngam
c***first executable statement  clbeta
      if (real(a) .le. 0.0 .or. real(b) .le. 0.0) call xermsg ('slatec',
     +   'clbeta', 'real part of both arguments must be gt 0', 1, 2)
c
      clbeta = clngam(a) + clngam(b) - clngam(a+b)
c
      return
      end
