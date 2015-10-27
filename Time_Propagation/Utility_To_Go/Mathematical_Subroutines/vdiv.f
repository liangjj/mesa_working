*deck @(#)vdiv.f	5.1  11/6/94
      subroutine vdiv(v,w,x,n)
c***begin prologue     vdiv
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, division
c***author             saxe, paul (lanl)
c***source             @(#)vdiv.f	5.1   11/6/94
c***purpose            vectorized vector division:  v=w/x .
c***description
c                      call vdiv(v,w,x,n)
c                        v        output vector of length (n).
c                        w        input vector of length (n).
c                        x        input vector of length (n).
c                        n        vector length.
c
c***references
c***routines called    (none)
c***end prologue       vdiv
      real*8 v(n),w(n),x(n)
      do 1 i=1,n
         v(i)=w(i)/x(i)
    1 continue
      return
      end
