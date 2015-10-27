*deck @(#)vsub.f	5.1  11/6/94
      subroutine vsub(v,w,x,n)
c***begin prologue     vsub
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, subtract
c***author             saxe, paul (lanl)
c***source             @(#)vsub.f	5.1   11/6/94
c***purpose            vectorized vector subtract:  v=w-x .
c***description
c                      call vsub(v,w,x,n)
c                        v        output vector of length n.
c                        w        input vector of length n.
c                        x        input vector of length n.
c                        n        vector lengths.
c
c***references
c***routines called    (none)
c***end prologue       vsub
      real*8 v(n),w(n),x(n)
      do 1 i=1,n
         v(i)=w(i)-x(i)
    1 continue
      return
      end
