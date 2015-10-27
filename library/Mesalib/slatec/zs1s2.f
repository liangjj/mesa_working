*deck zs1s2
      subroutine zs1s2 (zrr, zri, s1r, s1i, s2r, s2i, nz, ascle, alim,
     +   iuf)
c***begin prologue  zs1s2
c***subsidiary
c***purpose  subsidiary to zairy and zbesk
c***library   slatec
c***type      all (cs1s2-a, zs1s2-a)
c***author  amos, d. e., (snl)
c***description
c
c     zs1s2 tests for a possible underflow resulting from the
c     addition of the i and k functions in the analytic con-
c     tinuation formula where s1=k function and s2=i function.
c     on kode=1 the i and k functions are different orders of
c     magnitude, but for kode=2 they can be of the same order
c     of magnitude and the maximum must be at least one
c     precision above the underflow limit.
c
c***see also  zairy, zbesk
c***routines called  zabs, zexp, zlog
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c   930122  added zexp and zlog to external statement.  (rwc)
c***end prologue  zs1s2
c     complex czero,c1,s1,s1d,s2,zr
      double precision aa, alim, aln, ascle, as1, as2, c1i, c1r, s1di,
     * s1dr, s1i, s1r, s2i, s2r, zeroi, zeror, zri, zrr, zabs
      integer iuf, idum, nz
      external zabs, zexp, zlog
      data zeror,zeroi  / 0.0d0 , 0.0d0 /
c***first executable statement  zs1s2
      nz = 0
      as1 = zabs(s1r,s1i)
      as2 = zabs(s2r,s2i)
      if (s1r.eq.0.0d0 .and. s1i.eq.0.0d0) go to 10
      if (as1.eq.0.0d0) go to 10
      aln = -zrr - zrr + log(as1)
      s1dr = s1r
      s1di = s1i
      s1r = zeror
      s1i = zeroi
      as1 = zeror
      if (aln.lt.(-alim)) go to 10
      call zlog(s1dr, s1di, c1r, c1i, idum)
      c1r = c1r - zrr - zrr
      c1i = c1i - zri - zri
      call zexp(c1r, c1i, s1r, s1i)
      as1 = zabs(s1r,s1i)
      iuf = iuf + 1
   10 continue
      aa = max(as1,as2)
      if (aa.gt.ascle) return
      s1r = zeror
      s1i = zeroi
      s2r = zeror
      s2i = zeroi
      nz = 1
      iuf = 0
      return
      end
