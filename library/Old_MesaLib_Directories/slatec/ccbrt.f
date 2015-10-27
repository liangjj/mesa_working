*deck ccbrt
      complex function ccbrt (z)
c***begin prologue  ccbrt
c***purpose  compute the cube root.
c***library   slatec (fnlib)
c***category  c2
c***type      complex (cbrt-s, dcbrt-d, ccbrt-c)
c***keywords  cube root, elementary functions, fnlib, roots
c***author  fullerton, w., (lanl)
c***description
c
c ccbrt(z) calculates the complex cube root of z.  the principal root
c for which -pi .lt. arg(z) .le. +pi is returned.
c
c***references  (none)
c***routines called  carg, cbrt
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  ccbrt
      complex z
c***first executable statement  ccbrt
      theta = carg(z) / 3.0
      r = cbrt (abs(z))
c
      ccbrt = cmplx (r*cos(theta), r*sin(theta))
c
      return
      end
