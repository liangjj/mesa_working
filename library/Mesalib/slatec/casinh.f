*deck casinh
      complex function casinh (z)
c***begin prologue  casinh
c***purpose  compute the arc hyperbolic sine.
c***library   slatec (fnlib)
c***category  c4c
c***type      complex (asinh-s, dasinh-d, casinh-c)
c***keywords  arc hyperbolic sine, asinh, elementary functions, fnlib,
c             inverse hyperbolic sine
c***author  fullerton, w., (lanl)
c***description
c
c casinh(z) calculates the complex arc hyperbolic sine of z.
c
c***references  (none)
c***routines called  casin
c***revision history  (yymmdd)
c   770401  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  casinh
      complex z, ci, casin
      save ci
      data ci /(0.,1.)/
c***first executable statement  casinh
      casinh = -ci*casin (ci*z)
c
      return
      end
