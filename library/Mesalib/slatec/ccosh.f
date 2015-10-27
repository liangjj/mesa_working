*deck ccosh
      complex function ccosh (z)
c***begin prologue  ccosh
c***purpose  compute the complex hyperbolic cosine.
c***library   slatec (fnlib)
c***category  c4c
c***type      complex (ccosh-c)
c***keywords  elementary functions, fnlib, hyperbolic cosine
c***author  fullerton, w., (lanl)
c***description
c
c ccosh(z) calculates the complex hyperbolic cosine of z.
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  ccosh
      complex z, ci
      save ci
      data ci /(0.,1.)/
c***first executable statement  ccosh
      ccosh = cos (ci*z)
c
      return
      end
