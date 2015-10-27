*deck catanh
      complex function catanh (z)
c***begin prologue  catanh
c***purpose  compute the arc hyperbolic tangent.
c***library   slatec (fnlib)
c***category  c4c
c***type      complex (atanh-s, datanh-d, catanh-c)
c***keywords  arc hyperbolic tangent, atanh, elementary functions,
c             fnlib, inverse hyperbolic tangent
c***author  fullerton, w., (lanl)
c***description
c
c catanh(z) calculates the complex arc hyperbolic tangent of z.
c
c***references  (none)
c***routines called  catan
c***revision history  (yymmdd)
c   770401  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  catanh
      complex z, ci, catan
      save ci
      data ci /(0.,1.)/
c***first executable statement  catanh
      catanh = -ci*catan(ci*z)
c
      return
      end
