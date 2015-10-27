*deck cosdg
      function cosdg (x)
c***begin prologue  cosdg
c***purpose  compute the cosine of an argument in degrees.
c***library   slatec (fnlib)
c***category  c4a
c***type      single precision (cosdg-s, dcosdg-d)
c***keywords  cosine, degrees, elementary functions, fnlib,
c             trigonometric
c***author  fullerton, w., (lanl)
c***description
c
c cosdg(x) evaluates the cosine for real x in degrees.
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   770601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  cosdg
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
      save raddeg
      data raddeg / .017453292519943296e0 /
c
c***first executable statement  cosdg
      cosdg = cos (raddeg*x)
c
      if (mod(x,90.).ne.0.) return
      n = abs(x)/90.0 + 0.5
      n = mod (n, 2)
      if (n.eq.0) cosdg = sign (1.0, cosdg)
      if (n.eq.1) cosdg = 0.0
c
      return
      end
