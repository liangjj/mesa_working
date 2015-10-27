*deck wnlt3
      subroutine wnlt3 (i, imax, m, mdw, ipivot, h, w)
c***begin prologue  wnlt3
c***subsidiary
c***purpose  subsidiary to wnlit
c***library   slatec
c***type      single precision (wnlt3-s, dwnlt3-d)
c***author  hanson, r. j., (snla)
c           haskell, k. h., (snla)
c***description
c
c     perform column interchange.
c     exchange elements of permuted index vector and perform column
c     interchanges.
c
c***see also  wnlit
c***routines called  sswap
c***revision history  (yymmdd)
c   790701  date written
c   890620  code extracted from wnlt and made a subroutine.  (rwc))
c***end prologue  wnlt3
      integer i, imax, ipivot(*), m, mdw
      real             h(*), w(mdw,*)
c
      external sswap
c
      real             t
      integer itemp
c
c***first executable statement  wnlt3
      if (imax.ne.i) then
         itemp        = ipivot(i)
         ipivot(i)    = ipivot(imax)
         ipivot(imax) = itemp
c
         call sswap(m, w(1,imax), 1, w(1,i), 1)
c
         t       = h(imax)
         h(imax) = h(i)
         h(i)    = t
      endif
      return
      end
