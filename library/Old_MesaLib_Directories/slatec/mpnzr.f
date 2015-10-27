*deck mpnzr
      subroutine mpnzr (rs, re, z, trunc)
c***begin prologue  mpnzr
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpnzr-a)
c***author  (unknown)
c***description
c
c  modified for use with blas.  blank common changed to named common.
c  assumes long (i.e. (t+4)-digit) fraction in r, sign = rs, exponent
c  = re.  normalizes, and returns 'mp' result in z. integer arguments
c  rs and re are not preserved. r*-rounding is used if trunc.eq.0
c
c  the argument z(*) and the variable r in common are integer arrays
c  of size 30.  see the comments in the routine mpblas for the reason
c  for this choice.
c
c***see also  dqdota, dqdoti, mpblas
c***routines called  mperr, mpovfl, mpunfl
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   930124  increased array size in mpcon for sun -r8.  (rwc)
c***end prologue  mpnzr
      common /mpcom/ b, t, m, lun, mxr, r(30)
      integer b, t, r, z(*), re, rs, trunc, b2
c***first executable statement  mpnzr
      i2 = t + 4
      if (rs.ne.0) go to 20
c store zero in z
   10 z(1) = 0
      return
c check that sign = +-1
   20 if (abs(rs).le.1) go to 40
      write (lun, 30)
   30 format (' *** sign not 0, +1 or -1 in call to mpnzr,',
     1        ' possible overwriting problem ***')
      call mperr
      go to 10
c look for first nonzero digit
   40 do 50 i = 1, i2
      is = i - 1
      if (r(i).gt.0) go to 60
   50 continue
c fraction zero
      go to 10
   60 if (is.eq.0) go to 90
c normalize
      re = re - is
      i2m = i2 - is
      do 70 j = 1, i2m
      k = j + is
   70 r(j) = r(k)
      i2p = i2m + 1
      do 80 j = i2p, i2
   80 r(j) = 0
c check to see if truncation is desired
   90 if (trunc.ne.0) go to 150
c see if rounding necessary
c treat even and odd bases differently
      b2 = b/2
      if ((2*b2).ne.b) go to 130
c b even.  round if r(t+1).ge.b2 unless r(t) odd and all zeros
c after r(t+2).
      if (r(t+1) - b2) 150, 100, 110
  100 if (mod(r(t),2).eq.0) go to 110
      if ((r(t+2)+r(t+3)+r(t+4)).eq.0) go to 150
c round
  110 do 120 j = 1, t
      i = t + 1 - j
      r(i) = r(i) + 1
      if (r(i).lt.b) go to 150
  120 r(i) = 0
c exceptional case, rounded up to .10000...
      re = re + 1
      r(1) = 1
      go to 150
c odd base, round if r(t+1)... .gt. 1/2
  130 do 140 i = 1, 4
      it = t + i
      if (r(it) - b2) 150, 140, 110
  140 continue
c check for overflow
  150 if (re.le.m) go to 170
      write (lun, 160)
  160 format (' *** overflow occurred in mpnzr ***')
      call mpovfl (z)
      return
c check for underflow
  170 if (re.lt.(-m)) go to 190
c store result in z
      z(1) = rs
      z(2) = re
      do 180 i = 1, t
  180 z(i+2) = r(i)
      return
c underflow here
  190 call mpunfl (z)
      return
      end
