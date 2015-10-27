*deck @(#)vmove.f	5.1  11/6/94
      subroutine vmove(v,w,n)
c***begin prologue     vmove
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, move, transfer, copy
c***author             saxe, paul (lanl)
c***source             @(#)vmove.f	5.1   11/6/94
c***purpose            vectorized copy: v=w .
c***description
c                      call vmove(v,w,n)
c                        v        output vector of length n.
c                        w        input vector of length n.
c                        n        vector length.
c
c***references
c***routines called    (none)
c***end prologue       vmove
      real*8 v(n),w(n)
c
c      do 1 i=1,n
c         v(i)=w(i)
c    1 continue
      call scopy(n,w,1,v,1)
c
      return
      end
