*deck cuchk
      subroutine cuchk (y, nz, ascle, tol)
c***begin prologue  cuchk
c***subsidiary
c***purpose  subsidiary to seri, cuoik, cunk1, cunk2, cuni1, cuni2 and
c            ckscl
c***library   slatec
c***type      all (cuchk-a, zuchk-a)
c***author  amos, d. e., (snl)
c***description
c
c      y enters as a scaled quantity whose magnitude is greater than
c      exp(-alim)=ascle=1.0e+3*r1mach(1)/tol. the test is made to see
c      if the magnitude of the real or imaginary part would under flow
c      when y is scaled (by tol) to its proper value. y is accepted
c      if the underflow is at least one precision below the magnitude
c      of the largest component; otherwise the phase angle does not have
c      absolute accuracy and an underflow is assumed.
c
c***see also  ckscl, cuni1, cuni2, cunk1, cunk2, cuoik, seri
c***routines called  (none)
c***revision history  (yymmdd)
c   ??????  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cuchk
c
      complex y
      real ascle, ss, st, tol, yr, yi
      integer nz
c***first executable statement  cuchk
      nz = 0
      yr = real(y)
      yi = aimag(y)
      yr = abs(yr)
      yi = abs(yi)
      st = min(yr,yi)
      if (st.gt.ascle) return
      ss = max(yr,yi)
      st=st/tol
      if (ss.lt.st) nz = 1
      return
      end
