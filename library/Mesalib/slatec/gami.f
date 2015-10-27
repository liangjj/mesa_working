*deck gami
      function gami (a, x)
c***begin prologue  gami
c***purpose  evaluate the incomplete gamma function.
c***library   slatec (fnlib)
c***category  c7e
c***type      single precision (gami-s, dgami-d)
c***keywords  fnlib, incomplete gamma function, special functions
c***author  fullerton, w., (lanl)
c***description
c
c evaluate the incomplete gamma function defined by
c
c gami = integral from t = 0 to x of exp(-t) * t**(a-1.0) .
c
c gami is evaluated for positive values of a and non-negative values
c of x.  a slight deterioration of 2 or 3 digits accuracy will occur
c when gami is very large or very small, because logarithmic variables
c are used.  gami, a, and x are single precision.
c
c***references  (none)
c***routines called  alngam, gamit, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  gami
c***first executable statement  gami
      if (a .le. 0.0) call xermsg ('slatec', 'gami',
     +   'a must be gt zero', 1, 2)
      if (x .lt. 0.0) call xermsg ('slatec', 'gami',
     +   'x must be ge zero', 2, 2)
c
      gami = 0.0
      if (x.eq.0.0) return
c
c the only error possible in the expression below is a fatal overflow.
      factor = exp (alngam(a) + a*log(x) )
c
      gami = factor * gamit(a, x)
c
      return
      end
