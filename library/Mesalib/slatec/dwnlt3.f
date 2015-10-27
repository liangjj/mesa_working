*deck dwnlt3
      subroutine dwnlt3 (i, imax, m, mdw, ipivot, h, w)
c***begin prologue  dwnlt3
c***subsidiary
c***purpose  subsidiary to wnlit
c***library   slatec
c***type      double precision (wnlt3-s, dwnlt3-d)
c***author  hanson, r. j., (snla)
c           haskell, k. h., (snla)
c***description
c
c     perform column interchange.
c     exchange elements of permuted index vector and perform column
c     interchanges.
c
c***see also  dwnlit
c***routines called  dswap
c***revision history  (yymmdd)
c   790701  date written
c   890620  code extracted from wnlit and made a subroutine.  (rwc))
c   900604  dp version created from sp version.  (rwc)
c***end prologue  dwnlt3
      integer i, imax, ipivot(*), m, mdw
      double precision h(*), w(mdw,*)
c
      external dswap
c
      double precision t
      integer itemp
c
c***first executable statement  dwnlt3
      if (imax.ne.i) then
         itemp        = ipivot(i)
         ipivot(i)    = ipivot(imax)
         ipivot(imax) = itemp
c
         call dswap(m, w(1,imax), 1, w(1,i), 1)
c
         t       = h(imax)
         h(imax) = h(i)
         h(i)    = t
      endif
      return
      end
