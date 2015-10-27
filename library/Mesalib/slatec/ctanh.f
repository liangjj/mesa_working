*deck ctanh
      complex function ctanh (z)
c***begin prologue  ctanh
c***purpose  compute the complex hyperbolic tangent.
c***library   slatec (fnlib)
c***category  c4c
c***type      complex (ctanh-c)
c***keywords  elementary functions, fnlib, hyperbolic tangent
c***author  fullerton, w., (lanl)
c***description
c
c ctanh(z) calculates the complex hyperbolic tangent of complex
c argument z.  z is in units of radians.
c
c***references  (none)
c***routines called  ctan
c***revision history  (yymmdd)
c   770401  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  ctanh
      complex z, ci, ctan
      save ci
      data ci /(0.,1.)/
c***first executable statement  ctanh
      ctanh = -ci*ctan(ci*z)
c
      return
      end
