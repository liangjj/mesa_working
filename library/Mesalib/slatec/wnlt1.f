*deck wnlt1
      subroutine wnlt1 (i, lend, mend, ir, mdw, recalc, imax, hbar, h,
     +   scale, w)
c***begin prologue  wnlt1
c***subsidiary
c***purpose  subsidiary to wnlit
c***library   slatec
c***type      single precision (wnlt1-s, dwnlt1-d)
c***author  hanson, r. j., (snla)
c           haskell, k. h., (snla)
c***description
c
c     to update the column sum of squares and find the pivot column.
c     the column sum of squares vector will be updated at each step.
c     when numerically necessary, these values will be recomputed.
c
c***see also  wnlit
c***routines called  isamax
c***revision history  (yymmdd)
c   790701  date written
c   890620  code extracted from wnlit and made a subroutine.  (rwc))
c***end prologue  wnlt1
      integer i, imax, ir, lend, mdw, mend
      real             h(*), hbar, scale(*), w(mdw,*)
      logical recalc
c
      external isamax
      integer isamax
c
      integer j, k
c
c***first executable statement  wnlt1
      if (ir.ne.1 .and. (.not.recalc)) then
c
c        update column ss=sum of squares.
c
         do 10 j=i,lend
            h(j) = h(j) - scale(ir-1)*w(ir-1,j)**2
   10    continue
c
c        test for numerical accuracy.
c
         imax = isamax(lend-i+1, h(i), 1) + i - 1
         recalc = (hbar+1.e-3*h(imax)) .eq. hbar
      endif
c
c     if required, recalculate column ss, using rows ir through mend.
c
      if (recalc) then
         do 30 j=i,lend
            h(j) = 0.e0
            do 20 k=ir,mend
               h(j) = h(j) + scale(k)*w(k,j)**2
   20       continue
   30    continue
c
c        find column with largest ss.
c
         imax = isamax(lend-i+1, h(i), 1) + i - 1
         hbar = h(imax)
      endif
      return
      end
