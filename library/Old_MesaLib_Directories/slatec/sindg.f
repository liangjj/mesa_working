*deck sindg
      function sindg (x)
c***begin prologue  sindg
c***purpose  compute the sine of an argument in degrees.
c***library   slatec (fnlib)
c***category  c4a
c***type      single precision (sindg-s, dsindg-d)
c***keywords  degrees, elementary functions, fnlib, sine, trigonometric
c***author  fullerton, w., (lanl)
c***description
c
c sindg(x) evaluates the single precision sine of x where
c x is in degrees.
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  sindg
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
      save raddeg
      data raddeg / .017453292519943296e0 /
c
c***first executable statement  sindg
      sindg = sin (raddeg*x)
c
      if (mod(x,90.).ne.0.) return
      n = abs(x)/90.0 + 0.5
      n = mod (n, 2)
      if (n.eq.0) sindg = 0.
      if (n.eq.1) sindg = sign (1.0, sindg)
c
      return
      end
