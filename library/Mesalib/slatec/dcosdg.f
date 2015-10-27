*deck dcosdg
      double precision function dcosdg (x)
c***begin prologue  dcosdg
c***purpose  compute the cosine of an argument in degrees.
c***library   slatec (fnlib)
c***category  c4a
c***type      double precision (cosdg-s, dcosdg-d)
c***keywords  cosine, degrees, elementary functions, fnlib,
c             trigonometric
c***author  fullerton, w., (lanl)
c***description
c
c dcosdg(x) calculates the double precision trigonometric cosine
c for double precision argument x in units of degrees.
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  dcosdg
      double precision x, raddeg
      save raddeg
      data raddeg / 0.0174532925 1994329576 9236907684 886 d0 /
c***first executable statement  dcosdg
      dcosdg = cos (raddeg*x)
c
      if (mod(x,90.d0).ne.0.d0) return
      n = abs(x)/90.d0 + 0.5d0
      n = mod (n, 2)
      if (n.eq.0) dcosdg = sign (1.0d0, dcosdg)
      if (n.eq.1) dcosdg = 0.0d0
c
      return
      end
