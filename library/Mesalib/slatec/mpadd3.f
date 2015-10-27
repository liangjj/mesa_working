*deck mpadd3
      subroutine mpadd3 (x, y, s, med, re)
c***begin prologue  mpadd3
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpadd3-a)
c***author  (unknown)
c***description
c
c   called by mpadd2; does inner loops of addition
c
c   the arguments x(*) and y(*) and the variable r in common are all
c   integer arrays of size 30.  see the comments in the routine mpblas
c   for the reason for this choice.
c
c***see also  dqdota, dqdoti, mpblas
c***routines called  (none)
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   ??????  modified for use with blas.  blank common changed to named
c           common.  r given dimension 12.
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   930124  increased array size in mpcon for sun -r8.  (rwc)
c***end prologue  mpadd3
      common /mpcom/ b, t, m, lun, mxr, r(30)
      integer b, t, r, x(*), y(*), s, re, c, ted
c***first executable statement  mpadd3
      ted = t + med
      i2 = t + 4
      i = i2
      c = 0
c clear guard digits to right of x digits
   10 if (i.le.ted) go to 20
      r(i) = 0
      i = i - 1
      go to 10
   20 if (s.lt.0) go to 130
c here do addition, exponent(y) .ge. exponent(x)
      if (i.lt.t) go to 40
   30 j = i - med
      r(i) = x(j+2)
      i = i - 1
      if (i.gt.t) go to 30
   40 if (i.le.med) go to 60
      j = i - med
      c = y(i+2) + x(j+2) + c
      if (c.lt.b) go to 50
c carry generated here
      r(i) = c - b
      c = 1
      i = i - 1
      go to 40
c no carry generated here
   50 r(i) = c
      c = 0
      i = i - 1
      go to 40
   60 if (i.le.0) go to 90
      c = y(i+2) + c
      if (c.lt.b) go to 70
      r(i) = 0
      c = 1
      i = i - 1
      go to 60
   70 r(i) = c
      i = i - 1
c no carry possible here
   80 if (i.le.0) return
      r(i) = y(i+2)
      i = i - 1
      go to 80
   90 if (c.eq.0) return
c must shift right here as carry off end
      i2p = i2 + 1
      do 100 j = 2, i2
      i = i2p - j
  100 r(i+1) = r(i)
      r(1) = 1
      re = re + 1
      return
c here do subtraction, abs(y) .gt. abs(x)
  110 j = i - med
      r(i) = c - x(j+2)
      c = 0
      if (r(i).ge.0) go to 120
c borrow generated here
      c = -1
      r(i) = r(i) + b
  120 i = i - 1
  130 if (i.gt.t) go to 110
  140 if (i.le.med) go to 160
      j = i - med
      c = y(i+2) + c - x(j+2)
      if (c.ge.0) go to 150
c borrow generated here
      r(i) = c + b
      c = -1
      i = i - 1
      go to 140
c no borrow generated here
  150 r(i) = c
      c = 0
      i = i - 1
      go to 140
  160 if (i.le.0) return
      c = y(i+2) + c
      if (c.ge.0) go to 70
      r(i) = c + b
      c = -1
      i = i - 1
      go to 160
      end
