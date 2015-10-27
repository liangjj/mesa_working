*deck casin
      complex function casin (zinp)
c***begin prologue  casin
c***purpose  compute the complex arc sine.
c***library   slatec (fnlib)
c***category  c4a
c***type      complex (casin-c)
c***keywords  arc sine, elementary functions, fnlib, trigonometric
c***author  fullerton, w., (lanl)
c***description
c
c casin(zinp) calculates the complex trigonometric arc sine of zinp.
c the result is in units of radians, and the real part is in the first
c or fourth quadrant.
c
c***references  (none)
c***routines called  r1mach
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  casin
      complex zinp, z, z2, sqzp1, ci
      logical first
      save pi2, pi, ci, nterms, rmin, first
      data pi2 /1.5707963267 9489661923e0/
      data pi /3.1415926535 8979324e0/
      data ci /(0.,1.)/
      data first /.true./
c***first executable statement  casin
      if (first) then
c nterms = log(eps)/log(rmax)  where rmax = 0.1
         nterms = -0.4343*log(r1mach(3))
         rmin = sqrt (6.0*r1mach(3))
      endif
      first = .false.
c
      z = zinp
      r = abs (z)
      if (r.gt.0.1) go to 30
c
      casin = z
      if (r.lt.rmin) return
c
      casin = (0.0, 0.0)
      z2 = z*z
      do 20 i=1,nterms
        twoi = 2*(nterms-i) + 1
        casin = 1.0/twoi + twoi*casin*z2/(twoi+1.0)
 20   continue
      casin = z*casin
      return
c
 30   if (real(zinp).lt.0.0) z = -zinp
c
      sqzp1 = sqrt (z+1.0)
      if (aimag(sqzp1).lt.0.) sqzp1 = -sqzp1
      casin = pi2 - ci * log (z + sqzp1*sqrt(z-1.0))
c
      if (real(casin).gt.pi2) casin = pi - casin
      if (real(casin).le.(-pi2)) casin = -pi - casin
      if (real(zinp).lt.0.) casin = -casin
c
      return
      end
