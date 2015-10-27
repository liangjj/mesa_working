*deck @(#)vneg.f	5.1  11/6/94
      subroutine vneg(v,w,n)
c***begin prologue     vneg
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, negate
c***author             saxe, paul (lanl)
c***source             @(#)vneg.f	5.1   11/6/94
c***purpose            vectorized vector negate: v=-w .
c***description
c                      call vneg(v,w,n)
c                        v        output vector.
c                        w        input vector.
c                        n        vector length.
c
c***references
c***routines called    (none)
c***end prologue       vneg
      real*8 v(n),w(n)
      do 1 i=1,n
         v(i)=-w(i)
    1 continue
      return
      end
