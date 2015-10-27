*deck bdiff
      subroutine bdiff (l, v)
c***begin prologue  bdiff
c***subsidiary
c***purpose  subsidiary to bskin
c***library   slatec
c***type      single precision (bdiff-s, dbdiff-d)
c***author  amos, d. e., (snla)
c***description
c
c     bdiff computes the sum of b(l,k)*v(k)*(-1)**k where b(l,k)
c     are the binomial coefficients.  truncated sums are computed by
c     setting last part of the v vector to zero. on return, the binomial
c     sum is in v(l).
c
c***see also  bskin
c***routines called  (none)
c***revision history  (yymmdd)
c   820601  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  bdiff
      integer i, j, k, l
      real v
      dimension v(*)
c***first executable statement  bdiff
      if (l.eq.1) return
      do 20 j=2,l
        k = l
        do 10 i=j,l
          v(k) = v(k-1) - v(k)
          k = k - 1
   10   continue
   20 continue
      return
      end
