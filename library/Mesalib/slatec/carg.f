*deck carg
      function carg (z)
c***begin prologue  carg
c***purpose  compute the argument of a complex number.
c***library   slatec (fnlib)
c***category  a4a
c***type      complex (carg-c)
c***keywords  argument of a complex number, elementary functions, fnlib
c***author  fullerton, w., (lanl)
c***description
c
c carg(z) calculates the argument of the complex number z.  note
c that carg returns a real result.  if z = x+iy, then carg is atan(y/x),
c except when both x and y are zero, in which case the result
c will be zero.
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   770401  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  carg
      complex z
c***first executable statement  carg
      carg = 0.0
      if (real(z).ne.0. .or. aimag(z).ne.0.) carg =
     1  atan2 (aimag(z), real(z))
c
      return
      end
