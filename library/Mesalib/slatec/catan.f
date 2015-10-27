*deck catan
      complex function catan (z)
c***begin prologue  catan
c***purpose  compute the complex arc tangent.
c***library   slatec (fnlib)
c***category  c4a
c***type      complex (catan-c)
c***keywords  arc tangent, elementary functions, fnlib, trigonometric
c***author  fullerton, w., (lanl)
c***description
c
c catan(z) calculates the complex trigonometric arc tangent of z.
c the result is in units of radians, and the real part is in the first
c or fourth quadrant.
c
c***references  (none)
c***routines called  r1mach, xermsg
c***revision history  (yymmdd)
c   770801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  catan
      complex z, z2
      logical first
      save pi2, nterms, sqeps, rmin, rmax, first
      data pi2 / 1.5707963267 9489661923e0 /
      data first /.true./
c***first executable statement  catan
      if (first) then
c nterms = log(eps)/log(rbnd) where rbnd = 0.1
         nterms = -0.4343*log(r1mach(3)) + 1.0
         sqeps = sqrt(r1mach(4))
         rmin = sqrt (3.0*r1mach(3))
         rmax = 1.0/r1mach(3)
      endif
      first = .false.
c
      r = abs(z)
      if (r.gt.0.1) go to 30
c
      catan = z
      if (r.lt.rmin) return
c
      catan = (0.0, 0.0)
      z2 = z*z
      do 20 i=1,nterms
        twoi = 2*(nterms-i) + 1
        catan = 1.0/twoi - z2*catan
 20   continue
      catan = z*catan
      return
c
 30   if (r.gt.rmax) go to 50
      x = real(z)
      y = aimag(z)
      r2 = r*r
      if (r2 .eq. 1.0 .and. x .eq. 0.0) call xermsg ('slatec', 'catan',
     +   'z is +i or -i', 2, 2)
      if (abs(r2-1.0).gt.sqeps) go to 40
      if (abs(cmplx(1.0, 0.0)+z*z) .lt. sqeps) call xermsg ('slatec',
     +   'catan', 'answer lt half precision, z**2 close to -1', 1, 1)
c
 40   xans = 0.5*atan2(2.0*x, 1.0-r2)
      yans = 0.25*log((r2+2.0*y+1.0)/(r2-2.0*y+1.0))
      catan = cmplx (xans, yans)
      return
c
 50   catan = cmplx (pi2, 0.)
      if (real(z).lt.0.0) catan = cmplx(-pi2,0.0)
      return
c
      end
