*deck csinh
      complex function csinh (z)
c***begin prologue  csinh
c***purpose  compute the complex hyperbolic sine.
c***library   slatec (fnlib)
c***category  c4c
c***type      complex (csinh-c)
c***keywords  elementary functions, fnlib, hyperbolic sine
c***author  fullerton, w., (lanl)
c***description
c
c csinh(z) calculates the complex hyperbolic sine of complex
c argument z.  z is in units of radians.
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  csinh
      complex z, ci
      save ci
      data ci /(0.,1.)/
c***first executable statement  csinh
      csinh = -ci*sin(ci*z)
c
      return
      end
