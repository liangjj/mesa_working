*deck dsindg
      double precision function dsindg (x)
c***begin prologue  dsindg
c***purpose  compute the sine of an argument in degrees.
c***library   slatec (fnlib)
c***category  c4a
c***type      double precision (sindg-s, dsindg-d)
c***keywords  degrees, elementary functions, fnlib, sine, trigonometric
c***author  fullerton, w., (lanl)
c***description
c
c dsindg(x) calculates the double precision sine for double
c precision argument x where x is in degrees.
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  dsindg
      double precision x, raddeg
      save raddeg
      data raddeg / 0.0174532925 1994329576 9236907684 886 d0 /
c***first executable statement  dsindg
      dsindg = sin (raddeg*x)
c
      if (mod(x,90.d0).ne.0.d0) return
      n = abs(x)/90.d0 + 0.5d0
      n = mod (n, 2)
      if (n.eq.0) dsindg = 0.d0
      if (n.eq.1) dsindg = sign (1.0d0, dsindg)
c
      return
      end
