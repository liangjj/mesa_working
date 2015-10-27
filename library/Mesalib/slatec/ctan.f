*deck ctan
      complex function ctan (z)
c***begin prologue  ctan
c***purpose  compute the complex tangent.
c***library   slatec (fnlib)
c***category  c4a
c***type      complex (ctan-c)
c***keywords  elementary functions, fnlib, tangent, trigonometric
c***author  fullerton, w., (lanl)
c***description
c
c ctan(z) calculates the complex trigonometric tangent of complex
c argument z.  z is in units of radians.
c
c***references  (none)
c***routines called  r1mach, xerclr, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  ctan
      complex z
      save sqeps
      data sqeps /0./
c***first executable statement  ctan
      if (sqeps.eq.0.) sqeps = sqrt (r1mach(4))
c
      x2 = 2.0*real(z)
      y2 = 2.0*aimag(z)
c
      sn2x = sin (x2)
      call xerclr
c
      den = cos(x2) + cosh(y2)
      if (den .eq. 0.) call xermsg ('slatec', 'ctan',
     +   'tan is singular for input z (x is pi/2 or 3*pi/2 and y is 0)',
     +   2, 2)
c
      if (abs(den).gt.max(abs(x2),1.)*sqeps) go to 10
      call xerclr
      call xermsg ('slatec', 'ctan',
     +   'answer lt half precision, abs(x) too big or x too near ' //
     +   'pi/2 or 3*pi/2', 1, 1)
c
 10   ctan = cmplx (sn2x/den, sinh(y2)/den)
c
      return
      end
