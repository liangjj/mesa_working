*deck mpdivi
      subroutine mpdivi (x, iy, z)
c***begin prologue  mpdivi
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpdivi-a)
c***author  (unknown)
c***description
c
c  divides 'mp' x by the single-precision integer iy giving 'mp' z.
c  this is much faster than division by an 'mp' number.
c
c  the arguments x(*) and z(*), and the variable r in common are all
c  integer arrays of size 30.  see the comments in the routine mpblas
c  for the reason for this choice.
c
c***see also  dqdota, dqdoti, mpblas
c***routines called  mpchk, mperr, mpnzr, mpstr, mpunfl
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   ??????  modified for use with blas.  blank common changed to named
c           common.  r given dimension 12.
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   930124  increased array size in mpcon for sun -r8.  (rwc)
c***end prologue  mpdivi
      common /mpcom/ b, t, m, lun, mxr, r(30)
      integer b, t, r, x(*), z(*), rs, re, r1, c, c2, b2
c***first executable statement  mpdivi
      rs = x(1)
      j = iy
      if (j) 30, 10, 40
   10 write (lun, 20)
   20 format (' *** attempted division by zero in call to mpdivi ***')
      go to 230
   30 j = -j
      rs = -rs
   40 re = x(2)
c check for zero dividend
      if (rs.eq.0) go to 120
c check for division by b
      if (j.ne.b) go to 50
      call mpstr (x, z)
      if (re.le.(-m)) go to 240
      z(1) = rs
      z(2) = re - 1
      return
c check for division by 1 or -1
   50 if (j.ne.1) go to 60
      call mpstr (x, z)
      z(1) = rs
      return
   60 c = 0
      i2 = t + 4
      i = 0
c if j*b not representable as an integer have to simulate
c long division.   assume at least 16-bit word.
      b2 = max(8*b,32767/b)
      if (j.ge.b2) go to 130
c look for first nonzero digit in quotient
   70 i = i + 1
      c = b*c
      if (i.le.t) c = c + x(i+2)
      r1 = c/j
      if (r1) 210, 70, 80
c adjust exponent and get t+4 digits in quotient
   80 re = re + 1 - i
      r(1) = r1
      c = b*(c - j*r1)
      kh = 2
      if (i.ge.t) go to 100
      kh = 1 + t - i
      do 90 k = 2, kh
      i = i + 1
      c = c + x(i+2)
      r(k) = c/j
   90 c = b*(c - j*r(k))
      if (c.lt.0) go to 210
      kh = kh + 1
  100 do 110 k = kh, i2
      r(k) = c/j
  110 c = b*(c - j*r(k))
      if (c.lt.0) go to 210
c normalize and round result
  120 call mpnzr (rs, re, z, 0)
      return
c here need simulated double-precision division
  130 c2 = 0
      j1 = j/b
      j2 = j - j1*b
      j11 = j1 + 1
c look for first nonzero digit
  140 i = i + 1
      c = b*c + c2
      c2 = 0
      if (i.le.t) c2 = x(i+2)
      if (c-j1) 140, 150, 160
  150 if (c2.lt.j2) go to 140
c compute t+4 quotient digits
  160 re = re + 1 - i
      k = 1
      go to 180
c main loop for large abs(iy) case
  170 k = k + 1
      if (k.gt.i2) go to 120
      i = i + 1
c get approximate quotient first
  180 ir = c/j11
c now reduce so overflow does not occur
      iq = c - ir*j1
      if (iq.lt.b2) go to 190
c here iq*b would possibly overflow so increase ir
      ir = ir + 1
      iq = iq - j1
  190 iq = iq*b - ir*j2
      if (iq.ge.0) go to 200
c here iq negative so ir was too large
      ir = ir - 1
      iq = iq + j
  200 if (i.le.t) iq = iq + x(i+2)
      iqj = iq/j
c r(k) = quotient, c = remainder
      r(k) = iqj + ir
      c = iq - j*iqj
      if (c.ge.0) go to 170
c carry negative so overflow must have occurred
  210 call mpchk (1, 4)
      write (lun, 220)
  220 format (' *** integer overflow in mpdivi, b too large ***')
  230 call mperr
      z(1) = 0
      return
c underflow here
  240 call mpunfl(z)
      return
      end
