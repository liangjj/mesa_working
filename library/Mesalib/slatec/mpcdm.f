*deck mpcdm
      subroutine mpcdm (dx, z)
c***begin prologue  mpcdm
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpcdm-a)
c***author  (unknown)
c***description
c
c converts double-precision number dx to multiple-precision z.
c some numbers will not convert exactly on machines with base
c other than two, four or sixteen. this routine is not called
c by any other routine in 'mp', so may be omitted if double-
c precision is not available.
c
c the argument z(*) and the variable r in common are both integer
c arrays of size 30.  see the comments in the routine mpblas for the
c for the reason for this choice.
c
c***see also  dqdota, dqdoti, mpblas
c***routines called  mpchk, mpdivi, mpmuli, mpnzr
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   ??????  modified for use with blas.  blank common changed to named
c           common.  r given dimension 12.
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   930124  increased array size in mpcon for sun -r8.  (rwc)
c***end prologue  mpcdm
      double precision db, dj, dx
      common /mpcom/ b, t, m, lun, mxr, r(30)
      integer b, t, r, z(*), rs, re, tp
c***first executable statement  mpcdm
      call mpchk (1, 4)
      i2 = t + 4
c check sign
      if (dx) 20, 10, 30
c if dx = 0d0 return 0
   10 z(1) = 0
      return
c dx .lt. 0d0
   20 rs = -1
      dj = -dx
      go to 40
c dx .gt. 0d0
   30 rs = 1
      dj = dx
   40 ie = 0
   50 if (dj.lt.1d0) go to 60
c increase ie and divide dj by 16.
      ie = ie + 1
      dj = 0.0625d0*dj
      go to 50
   60 if (dj.ge.0.0625d0) go to 70
      ie = ie - 1
      dj = 16d0*dj
      go to 60
c now dj is dy divided by suitable power of 16
c set exponent to 0
   70 re = 0
      db = dble(b)
c conversion loop (assume double-precision ops. exact)
      do 80 i = 1, i2
      dj = db*dj
      r(i) = int(dj)
   80 dj = dj - dble(r(i))
c normalize result
      call mpnzr (rs, re, z, 0)
      ib = max(7*b*b, 32767)/16
      tp = 1
c now multiply by 16**ie
      if (ie) 90, 130, 110
   90 k = -ie
      do 100 i = 1, k
      tp = 16*tp
      if ((tp.le.ib).and.(tp.ne.b).and.(i.lt.k)) go to 100
      call mpdivi (z, tp, z)
      tp = 1
  100 continue
      return
  110 do 120 i = 1, ie
      tp = 16*tp
      if ((tp.le.ib).and.(tp.ne.b).and.(i.lt.ie)) go to 120
      call mpmuli (z, tp, z)
      tp = 1
  120 continue
  130 return
      end
