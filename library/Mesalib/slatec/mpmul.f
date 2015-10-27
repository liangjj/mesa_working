*deck mpmul
      subroutine mpmul (x, y, z)
c***begin prologue  mpmul
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpmul-a)
c***author  (unknown)
c***description
c
c  multiplies x and y, returning result in z, for 'mp' x, y and z.
c  the simple o(t**2) algorithm is used, with four guard digits and
c  r*-rounding. advantage is taken of zero digits in x, but not in y.
c  asymptotically faster algorithms are known (see knuth, vol. 2),
c  but are difficult to implement in fortran in an efficient and
c  machine-independent manner. in comments to other 'mp' routines,
c  m(t) is the time to perform t-digit 'mp' multiplication. thus
c  m(t) = o(t**2) with the present version of mpmul, but
c  m(t) = o(t.log(t).log(log(t))) is theoretically possible.
c
c  the arguments x(*), y(*), and z(*), and the variable r in common are
c  all integer arrays of size 30.  see the comments in the routine
c  mpblas for the reason for this choice.
c
c***see also  dqdota, dqdoti, mpblas
c***routines called  mpchk, mperr, mpmlp, mpnzr
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   ??????  modified for use with blas.  blank common changed to named
c           common.  r given dimension 12.
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   930124  increased array size in mpcon for sun -r8.  (rwc)
c***end prologue  mpmul
      common /mpcom/ b, t, m, lun, mxr, r(30)
      integer b, t, r, x(*), y(*), z(*), rs, re, xi, c, ri
c***first executable statement  mpmul
      call mpchk (1, 4)
      i2 = t + 4
      i2p = i2 + 1
c form sign of product
      rs = x(1)*y(1)
      if (rs.ne.0) go to 10
c set result to zero
      z(1) = 0
      return
c form exponent of product
   10 re = x(2) + y(2)
c clear accumulator
      do 20 i = 1, i2
   20 r(i) = 0
c perform multiplication
      c = 8
      do 40 i = 1, t
      xi = x(i+2)
c for speed, put the number with many zeros first
      if (xi.eq.0) go to 40
      call mpmlp (r(i+1), y(3), xi, min (t, i2 - i))
      c = c - 1
      if (c.gt.0) go to 40
c check for legal base b digit
      if ((xi.lt.0).or.(xi.ge.b)) go to 90
c propagate carries at end and every eighth time,
c faster than doing it every time.
      do 30 j = 1, i2
      j1 = i2p - j
      ri = r(j1) + c
      if (ri.lt.0) go to 70
      c = ri/b
   30 r(j1) = ri - b*c
      if (c.ne.0) go to 90
      c = 8
   40 continue
      if (c.eq.8) go to 60
      if ((xi.lt.0).or.(xi.ge.b)) go to 90
      c = 0
      do 50 j = 1, i2
      j1 = i2p - j
      ri = r(j1) + c
      if (ri.lt.0) go to 70
      c = ri/b
   50 r(j1) = ri - b*c
      if (c.ne.0) go to 90
c normalize and round result
   60 call mpnzr (rs, re, z, 0)
      return
   70 write (lun, 80)
   80 format (' *** integer overflow in mpmul, b too large ***')
      go to 110
   90 write (lun, 100)
  100 format (' *** illegal base b digit in call to mpmul,',
     1        ' possible overwriting problem ***')
  110 call mperr
      z(1) = 0
      return
      end
