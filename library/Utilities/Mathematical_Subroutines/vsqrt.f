*deck @(#)vsqrt.f	5.1  11/6/94
      subroutine vsqrt(v,w,n)
c***begin prologue     vsqrt
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, square root.
c***author             saxe, paul (lanl)
c***source             @(#)vsqrt.f	5.1   11/6/94
c***purpose            vectorized vector square root: v=sqrt(w) .
c***description
c                      call vsqrt(v,w,n)
c                        v        output vector dimensioned (n).
c                        w        input vector dimensioned (n).
c                        n        vector length.
c
c***references
c***routines called    (none)
c***end prologue       vsqrt
      real*8 v(n),w(n)
      do 1 i=1,n
         v(i)=sqrt(w(i))
    1 continue
      return
      end
