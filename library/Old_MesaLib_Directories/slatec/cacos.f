*deck cacos
      complex function cacos (z)
c***begin prologue  cacos
c***purpose  compute the complex arc cosine.
c***library   slatec (fnlib)
c***category  c4a
c***type      complex (cacos-c)
c***keywords  arc cosine, elementary functions, fnlib, trigonometric
c***author  fullerton, w., (lanl)
c***description
c
c cacos(z) calculates the complex trigonometric arc cosine of z.
c the result is in units of radians, and the real part is in the
c first or second quadrant.
c
c***references  (none)
c***routines called  casin
c***revision history  (yymmdd)
c   770401  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  cacos
      complex z, casin
      save pi2
      data pi2 /1.5707963267 9489661923e0/
c***first executable statement  cacos
      cacos = pi2 - casin (z)
c
      return
      end
