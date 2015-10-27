*deck mpmul2
      subroutine mpmul2 (x, iy, z, trunc)
c***begin prologue  mpmul2
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpmul2-a)
c***author  (unknown)
c***description
c
c  multiplies 'mp' x by single-precision integer iy giving 'mp' z.
c  multiplication by 1 may be used to normalize a number even if some
c  digits are greater than b-1. result is rounded if trunc.eq.0,
c  otherwise truncated.
c
c  the arguments x(*) and z(*), and the variable r in common are all
c  integer arrays of size 30.  see the comments in the routine mpblas
c  for the reason for this choice.
c
c***see also  dqdota, dqdoti, mpblas
c***routines called  mpchk, mperr, mpnzr, mpovfl, mpstr
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   ??????  modified for use with blas.  blank common changed to named
c           common.  r given dimension 12.
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   930124  increased array size in mpcon for sun -r8.  (rwc)
c***end prologue  mpmul2
      common /mpcom/ b, t, m, lun, mxr, r(30)
      integer b, t, r, x(*), z(*), trunc, re, rs
      integer c, c1, c2, ri, t1, t3, t4
c***first executable statement  mpmul2
      rs = x(1)
      if (rs.eq.0) go to 10
      j = iy
      if (j) 20, 10, 50
c result zero
   10 z(1) = 0
      return
   20 j = -j
      rs = -rs
c check for multiplication by b
      if (j.ne.b) go to 50
      if (x(2).lt.m) go to 40
      call mpchk (1, 4)
      write (lun, 30)
   30 format (' *** overflow occurred in mpmul2 ***')
      call mpovfl (z)
      return
   40 call mpstr (x, z)
      z(1) = rs
      z(2) = x(2) + 1
      return
c set exponent to exponent(x) + 4
   50 re = x(2) + 4
c form product in accumulator
      c = 0
      t1 = t + 1
      t3 = t + 3
      t4 = t + 4
c if j*b not representable as an integer we have to simulate
c double-precision multiplication.
      if (j.ge.max(8*b, 32767/b)) go to 110
      do 60 ij = 1, t
      i = t1 - ij
      ri = j*x(i+2) + c
      c = ri/b
   60 r(i+4) = ri - b*c
c check for integer overflow
      if (ri.lt.0) go to 130
c have to treat first four words of r separately
      do 70 ij = 1, 4
      i = 5 - ij
      ri = c
      c = ri/b
   70 r(i) = ri - b*c
      if (c.eq.0) go to 100
c have to shift right here as carry off end
   80 do 90 ij = 1, t3
      i = t4 - ij
   90 r(i+1) = r(i)
      ri = c
      c = ri/b
      r(1) = ri - b*c
      re = re + 1
      if (c) 130, 100, 80
c normalize and round or truncate result
  100 call mpnzr (rs, re, z, trunc)
      return
c here j is too large for single-precision multiplication
  110 j1 = j/b
      j2 = j - j1*b
c form product
      do 120 ij = 1, t4
      c1 = c/b
      c2 = c - b*c1
      i = t1 - ij
      ix = 0
      if (i.gt.0) ix = x(i+2)
      ri = j2*ix + c2
      is = ri/b
      c = j1*ix + c1 + is
  120 r(i+4) = ri - b*is
      if (c) 130, 100, 80
c can only get here if integer overflow occurred
  130 call mpchk (1, 4)
      write (lun, 140)
  140 format (' *** integer overflow in mpmul2, b too large ***')
      call mperr
      go to 10
      end
