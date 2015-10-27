*deck @(#)vadmul.f	5.1  11/6/94
      subroutine vadmul(v,y,w,x,n)
c***begin prologue     vadmul
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, multiply
c***author             schneider, barry (nsf)
c***source             @(#)vadmul.f	5.1   11/6/94
c***purpose            vectorized vector multiply:  v= y + w*x .
c***description
c                      call vadmul(v,y,w,x,n)
c                        v        output vector of length n.
c                        w        input vector of length n.
c                        x        input vector of length n.
c                        y        input vector of length n.
c                        n        length of vector.
c
c***references
c***routines called    (none)
c***end prologue       vadmul
      real*8 v(n),w(n),x(n),y(n)
      do 1 i=1,n
         v(i) = y(i) + w(i)*x(i)
    1 continue
      return
      end
