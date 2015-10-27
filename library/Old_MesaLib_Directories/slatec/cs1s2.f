*deck cs1s2
      subroutine cs1s2 (zr, s1, s2, nz, ascle, alim, iuf)
c***begin prologue  cs1s2
c***subsidiary
c***purpose  subsidiary to cairy and cbesk
c***library   slatec
c***type      all (cs1s2-a, zs1s2-a)
c***author  amos, d. e., (snl)
c***description
c
c     cs1s2 tests for a possible underflow resulting from the
c     addition of the i and k functions in the analytic con-
c     tinuation formula where s1=k function and s2=i function.
c     on kode=1 the i and k functions are different orders of
c     magnitude, but for kode=2 they can be of the same order
c     of magnitude and the maximum must be at least one
c     precision above the underflow limit.
c
c***see also  cairy, cbesk
c***routines called  (none)
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cs1s2
      complex czero, c1, s1, s1d, s2, zr
      real aa, alim, aln, ascle, as1, as2, xx
      integer iuf, nz
      data czero / (0.0e0,0.0e0) /
c***first executable statement  cs1s2
      nz = 0
      as1 = abs(s1)
      as2 = abs(s2)
      aa = real(s1)
      aln = aimag(s1)
      if (aa.eq.0.0e0 .and. aln.eq.0.0e0) go to 10
      if (as1.eq.0.0e0) go to 10
      xx = real(zr)
      aln = -xx - xx + alog(as1)
      s1d = s1
      s1 = czero
      as1 = 0.0e0
      if (aln.lt.(-alim)) go to 10
      c1 = clog(s1d) - zr - zr
      s1 = cexp(c1)
      as1 = abs(s1)
      iuf = iuf + 1
   10 continue
      aa = max(as1,as2)
      if (aa.gt.ascle) return
      s1 = czero
      s2 = czero
      nz = 1
      iuf = 0
      return
      end
