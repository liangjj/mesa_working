*deck catan
      function catan(z)
c***begin prologue  catan
c***date written   770801   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c4a
c***keywords  arc tangent,complex,elementary function
c***author  fullerton, w., (lanl)
c***purpose  computes the complex arc tangent.
c***description
c
c catan(z) calculates the complex trigonometric arc tangent of z.
c the result is in units of radians, and the real part is in the first
c or fourth quadrant.
c***references  (none)
c***routines called  r1mach,lnkerr
c***end prologue  catan
      implicit real*8 (a-h,o-z)
      complex*16 z, z2, catan, czero
      data pi2 / 1.57079632679489661923d0 /
      data nterms, sqeps, rmin, rmax / 0, 3*0.0d+00 /
c***first executable statement  catan
      zero=0.d0
      czero=dcmplx(0.d0,0.d0)
      if (nterms.ne.0) go to 10
c nterms = log(eps)/log(rbnd) where rbnd = 0.1
      nterms = -0.4343d+00*log(r1mach(3)) + 1.0d+00
      sqeps = sqrt(r1mach(4))
      rmin = sqrt (3.0d+00*r1mach(3))
      rmax = 1.0d+00/r1mach(3)
c
 10   r = abs(z)
      if (r.gt.0.1d+00) go to 30
c
      catan = z
      if (r.lt.rmin) return
c
      catan = czero
      z2 = z*z
      do 20 i=1,nterms
        twoi = 2*(nterms-i) + 1
        catan = 1.0d+00/twoi - z2*catan
 20   continue
      catan = z*catan
      return
c
 30   if (r.gt.rmax) go to 50
      x = real(z)
      y = imag(z)
      r2 = r*r
      if (r2.eq.1.0d+00 .and. x.eq.0.0d+00) then
           call lnkerr ( 'catan   z is +i or -i')
      endif
      if (abs(r2-1.0d0).gt.sqeps) go to 40
      if (abs(dcmplx(1.0d+00,0.0d+00)+z*z).lt.sqeps) then
          call lnkerr ( 'catan  answer lt half precision, z**2 close')
      endif
c
 40   xans = 0.5d+00*atan2(2.0d+00*x, 1.0d+00-r2)
      yans = 0.25d+00*log((r2+2.0d+00*y+1.0d+00)/
     1                     (r2-2.0d+00*y+1.0d+00))
      catan = dcmplx (xans, yans)
      return
c
 50   catan = dcmplx (pi2, zero)
      if (real(z).lt.0.0d+00) then
          catan = dcmplx(-pi2,zero)
      endif
      return
c
      end
