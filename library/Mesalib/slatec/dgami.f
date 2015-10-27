*deck dgami
      double precision function dgami (a, x)
c***begin prologue  dgami
c***purpose  evaluate the incomplete gamma function.
c***library   slatec (fnlib)
c***category  c7e
c***type      double precision (gami-s, dgami-d)
c***keywords  fnlib, incomplete gamma function, special functions
c***author  fullerton, w., (lanl)
c***description
c
c evaluate the incomplete gamma function defined by
c
c dgami = integral from t = 0 to x of exp(-t) * t**(a-1.0) .
c
c dgami is evaluated for positive values of a and non-negative values
c of x.  a slight deterioration of 2 or 3 digits accuracy will occur
c when dgami is very large or very small, because logarithmic variables
c are used.  the function and both arguments are double precision.
c
c***references  (none)
c***routines called  dgamit, dlngam, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dgami
      double precision a, x, factor, dlngam, dgamit
c***first executable statement  dgami
      if (a .le. 0.d0) call xermsg ('slatec', 'dgami',
     +   'a must be gt zero', 1, 2)
      if (x .lt. 0.d0) call xermsg ('slatec', 'dgami',
     +   'x must be ge zero', 2, 2)
c
      dgami = 0.d0
      if (x.eq.0.0d0) return
c
c the only error possible in the expression below is a fatal overflow.
      factor = exp (dlngam(a) + a*log(x))
c
      dgami = factor * dgamit (a, x)
c
      return
      end
