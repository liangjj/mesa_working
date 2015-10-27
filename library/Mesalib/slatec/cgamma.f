*deck cgamma
      complex function cgamma (z)
c***begin prologue  cgamma
c***purpose  compute the complete gamma function.
c***library   slatec (fnlib)
c***category  c7a
c***type      complex (gamma-s, dgamma-d, cgamma-c)
c***keywords  complete gamma function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c cgamma(z) calculates the complete gamma function for complex
c argument z.  this is a preliminary version that is portable
c but not accurate.
c
c***references  (none)
c***routines called  clngam
c***revision history  (yymmdd)
c   770701  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  cgamma
      complex z, clngam
c***first executable statement  cgamma
      cgamma = exp (clngam(z))
c
      return
      end
