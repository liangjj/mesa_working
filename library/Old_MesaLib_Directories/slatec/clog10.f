*deck clog10
      complex function clog10 (z)
c***begin prologue  clog10
c***purpose  compute the principal value of the complex base 10
c            logarithm.
c***library   slatec (fnlib)
c***category  c4b
c***type      complex (clog10-c)
c***keywords  base ten logarithm, elementary functions, fnlib
c***author  fullerton, w., (lanl)
c***description
c
c clog10(z) calculates the principal value of the complex common
c or base 10 logarithm of z for -pi .lt. arg(z) .le. +pi.
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  clog10
      complex z
      save aloge
      data aloge / 0.4342944819 0325182765e0 /
c***first executable statement  clog10
      clog10 = aloge * log(z)
c
      return
      end
