*deck zuchk
      subroutine zuchk (yr, yi, nz, ascle, tol)
c***begin prologue  zuchk
c***subsidiary
c***purpose  subsidiary to seri, zuoik, zunk1, zunk2, zuni1, zuni2 and
c            zkscl
c***library   slatec
c***type      all (cuchk-a, zuchk-a)
c***author  amos, d. e., (snl)
c***description
c
c      y enters as a scaled quantity whose magnitude is greater than
c      exp(-alim)=ascle=1.0e+3*d1mach(1)/tol. the test is made to see
c      if the magnitude of the real or imaginary part would underflow
c      when y is scaled (by tol) to its proper value. y is accepted
c      if the underflow is at least one precision below the magnitude
c      of the largest component; otherwise the phase angle does not have
c      absolute accuracy and an underflow is assumed.
c
c***see also  seri, zkscl, zuni1, zuni2, zunk1, zunk2, zuoik
c***routines called  (none)
c***revision history  (yymmdd)
c   ??????  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zuchk
c
c     complex y
      double precision ascle, ss, st, tol, wr, wi, yr, yi
      integer nz
c***first executable statement  zuchk
      nz = 0
      wr = abs(yr)
      wi = abs(yi)
      st = min(wr,wi)
      if (st.gt.ascle) return
      ss = max(wr,wi)
      st = st/tol
      if (ss.lt.st) nz = 1
      return
      end
