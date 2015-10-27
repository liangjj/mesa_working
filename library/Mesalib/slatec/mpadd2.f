*deck mpadd2
      subroutine mpadd2 (x, y, z, y1, trunc)
c***begin prologue  mpadd2
c***subsidiary
c***purpose  subsidiary to dqdota and dqdoti
c***library   slatec
c***type      all (mpadd2-a)
c***author  (unknown)
c***description
c
c  called by mpadd, mpsub etc.
c  x, y and z are mp numbers, y1 and trunc are integers.
c  to force call by reference rather than value/result, y1 is
c  declared as an array, but only y1(1) is ever used.
c  sets z = x + y1(1)*abs(y), where y1(1) = +- y(1).
c  if trunc .eq. 0, r*-rounding is used;  otherwise, truncation.
c  r*-rounding is defined in the kuki and cody reference.
c
c  the arguments x(*), y(*), and z(*) are all integer arrays of size
c  30.  see the comments in the routine mpblas for the reason for this
c  choice.
c
c***see also  dqdota, dqdoti, mpblas
c***references  h. kuki and w. j. cody, a statistical study of floating
c                 point number systems, communications of the acm 16, 4
c                 (april 1973), pp. 223-230.
c               r. p. brent, on the precision attainable with various
c                 floating-point number systems, ieee transactions on
c                 computers c-22, 6 (june 1973), pp. 601-607.
c               r. p. brent, a fortran multiple-precision arithmetic
c                 package, acm transactions on mathematical software 4,
c                 1 (march 1978), pp. 57-70.
c               r. p. brent, mp, a fortran multiple-precision arithmetic
c                 package, algorithm 524, acm transactions on mathema-
c                 tical software 4, 1 (march 1978), pp. 71-81.
c***routines called  mpadd3, mpchk, mperr, mpnzr, mpstr
c***common blocks    mpcom
c***revision history  (yymmdd)
c   791001  date written
c   ??????  modified for use with blas.  blank common changed to named
c           common.  r given dimension 12.
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c   920528  added a references section revised.  (wrb)
c   930124  increased array size in mpcon for sun -r8.  (rwc)
c***end prologue  mpadd2
      common /mpcom/ b, t, m, lun, mxr, r(30)
      integer b, t, r, x(*), y(*), z(*), y1(*), trunc
      integer s, ed, rs, re
c***first executable statement  mpadd2
      if (x(1).ne.0) go to 20
   10 call mpstr(y, z)
      z(1) = y1(1)
      return
   20 if (y1(1).ne.0) go to 40
   30 call mpstr (x, z)
      return
c compare signs
   40 s = x(1)*y1(1)
      if (abs(s).le.1) go to 60
      call mpchk (1, 4)
      write (lun, 50)
   50 format (' *** sign not 0, +1 or -1 in call to mpadd2,',
     1        ' possible overwriting problem ***')
      call mperr
      z(1) = 0
      return
c compare exponents
   60 ed = x(2) - y(2)
      med = abs(ed)
      if (ed) 90, 70, 120
c exponents equal so compare signs, then fractions if nec.
   70 if (s.gt.0) go to 100
      do 80 j = 1, t
      if (x(j+2) - y(j+2)) 100, 80, 130
   80 continue
c result is zero
      z(1) = 0
      return
c here exponent(y) .ge. exponent(x)
   90 if (med.gt.t) go to 10
  100 rs = y1(1)
      re = y(2)
      call mpadd3 (x, y, s, med, re)
c normalize, round or truncate, and return
  110 call mpnzr (rs, re, z, trunc)
      return
c abs(x) .gt. abs(y)
  120 if (med.gt.t) go to 30
  130 rs = x(1)
      re = x(2)
      call mpadd3 (y, x, s, med, re)
      go to 110
      end
