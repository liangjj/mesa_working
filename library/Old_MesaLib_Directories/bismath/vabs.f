*deck @(#)vabs.f	1.1  12/11/94
      subroutine vabs(v,w,n)
c***begin prologue     vabs
c***date written       941211  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, absolute
c***author             schneider, barry (nsf)
c***source             @(#)vabs.f
c***purpose            vectorized absolute value:  v=abs(w) .
c***description
c                      call vabs(v,w,n)
c                        v        output vector of length n.
c                        w        input vector of length n.
c                        n        vector lengths.
c
c***references
c***routines called    (none)
c***end prologue       vabs
      real*8 v(n),w(n)
      do 1 i=1,n
         v(i)=abs(w(i))
    1 continue
      return
      end
