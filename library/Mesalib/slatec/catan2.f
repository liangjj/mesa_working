*deck catan2
      complex function catan2 (csn, ccs)
c***begin prologue  catan2
c***purpose  compute the complex arc tangent in the proper quadrant.
c***library   slatec (fnlib)
c***category  c4a
c***type      complex (catan2-c)
c***keywords  arc tangent, elementary functions, fnlib, polar angel,
c             quadrant, trigonometric
c***author  fullerton, w., (lanl)
c***description
c
c catan2(csn,ccs) calculates the complex trigonometric arc
c tangent of the ratio csn/ccs and returns a result whose real
c part is in the correct quadrant (within a multiple of 2*pi).  the
c result is in units of radians and the real part is between -pi
c and +pi.
c
c***references  (none)
c***routines called  catan, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  catan2
      complex csn, ccs, catan
      save pi
      data pi / 3.1415926535 8979323846e0 /
c***first executable statement  catan2
      if (abs(ccs).eq.0.) go to 10
c
      catan2 = catan (csn/ccs)
      if (real(ccs).lt.0.) catan2 = catan2 + pi
      if (real(catan2).gt.pi) catan2 = catan2 - 2.0*pi
      return
c
 10   if (abs(csn) .eq. 0.) call xermsg ('slatec', 'catan2',
     +   'called with both arguments zero', 1, 2)
c
      catan2 = cmplx (sign(0.5*pi,real(csn)), 0.0)
c
      return
      end
