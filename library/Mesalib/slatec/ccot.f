*deck ccot
      complex function ccot (z)
c***begin prologue  ccot
c***purpose  compute the cotangent.
c***library   slatec (fnlib)
c***category  c4a
c***type      complex (cot-s, dcot-d, ccot-c)
c***keywords  cotangent, elementary functions, fnlib, trigonometric
c***author  fullerton, w., (lanl)
c***description
c
c ccot(z) calculates the complex trigonometric cotangent of z.
c
c***references  (none)
c***routines called  r1mach, xerclr, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  ccot
      complex z
      save sqeps
      data sqeps /0./
c***first executable statement  ccot
      if (sqeps.eq.0.) sqeps = sqrt (r1mach(4))
c
      x2 = 2.0*real(z)
      y2 = 2.0*aimag(z)
c
      sn2x = sin (x2)
      call xerclr
c
      den = cosh(y2) - cos(x2)
      if (den .eq. 0.) call xermsg ('slatec', 'ccot',
     +   'cot is singular for input z (x is 0 or pi and y is 0)', 2, 2)
c
      if (abs(den).gt.max(abs(x2),1.)*sqeps) go to 10
      call xerclr
      call xermsg ('slatec', 'ccot',
     +   'answer lt half precision, abs(x) too big or x too near ' //
     +   '0 or pi', 1, 1)
c
 10   ccot = cmplx (sn2x/den, -sinh(y2)/den)
c
      return
      end
