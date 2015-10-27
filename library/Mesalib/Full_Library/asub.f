*deck @(#)asub.f	1.1  11/30/90
      subroutine asub(n,a,b,c)
c***begin prologue     asub
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, subtract
c***author             martin, richard (lanl)
c***source             @(#)asub.f	1.1   11/30/90
c***purpose            vectorized vector subtract:  c=a-b.
c***description
c                      call asub(n,a,b,c)
c                        n        vector lengths.
c                        a        input vector of length n.
c                        b        input vector of length n.
c                        c        output vector of length n.
c
c***references
c***routines called    (none)
c***end prologue       asub
      implicit real*8(a-h,o-z)
      dimension c(1),a(1),b(1)
c
c
      if(n.lt.1) return
      do 1 i=1,n
    1     c(i)=a(i)-b(i)
      return
      end
