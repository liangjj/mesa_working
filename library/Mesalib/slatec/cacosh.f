*deck cacosh
      complex function cacosh (z)
c***begin prologue  cacosh
c***purpose  compute the arc hyperbolic cosine.
c***library   slatec (fnlib)
c***category  c4c
c***type      complex (acosh-s, dacosh-d, cacosh-c)
c***keywords  acosh, arc hyperbolic cosine, elementary functions, fnlib,
c             inverse hyperbolic cosine
c***author  fullerton, w., (lanl)
c***description
c
c cacosh(z) calculates the complex arc hyperbolic cosine of z.
c
c***references  (none)
c***routines called  cacos
c***revision history  (yymmdd)
c   770401  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  cacosh
      complex z, ci, cacos
      save ci
      data ci /(0.,1.)/
c***first executable statement  cacosh
      cacosh = ci*cacos(z)
c
      return
      end
