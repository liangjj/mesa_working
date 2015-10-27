      function ccot(z)
c***begin prologue  ccot
c***date written   770401   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c4a
c***keywords  complex,cotangent,elementary function
c***author  fullerton, w., (lanl)
c***purpose  computes the complex cotangent.
c***description
c
c ccot(z) calculates the comlex trigonometric cotangent of z.
c***references  (none)
c***routines called  r1mach,xerclr,xerror
c***end prologue  ccot
      implicit real*8 (a-h,o-z)
      complex*16 ccot, z
      data sqeps /0.d0/
c***first executable statement  ccot
      if (sqeps.eq.0.d0) sqeps = sqrt (r1mach(4))
c
      x2 = 2.d0*real(z)
      y2 = 2.d0*imag(z)
c
      sn2x = sin (x2)
      call xerclr
c
      den = cosh(y2) - cos(x2)
      if (den.eq.0.d0) call lnkerr (  'ccot    cot is singular for input
     1 z (x is 0 or pi and y is 0)')
c
      if (abs(den).gt.max(abs(x2),1.d0)*sqeps) go to 10
      call xerclr
      call lnkerr ( 'ccot    answer lt half precision, abs(x) too big or
     1 x too near 0 or pi')
c
 10   ccot = cmplx (sn2x/den, -sinh(y2)/den)
c
      return
      end
