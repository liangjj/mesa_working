*deck @(#)gvadd.f	5.1  11/6/94
      subroutine gvadd(v,nv,w,nw,x,nx,n)
c***begin prologue     gvadd
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector addition, variable stride.
c***author             saxe, paul (lanl)
c***source             @(#)gvadd.f	5.1   11/6/94
c***purpose            vectorized vector addition:  v=w+x .
c***description
c                      call gvadd(v,nv,w,nw,x,nx,n)
c                        v        input vector of length (nv*n).
c                        w        input vector of length (nw*n).
c                        x        input vector of length (nx*n).
c                        nv       stride on v.
c                        nw       stride on w.
c                        nx       stride on x.
c                        n        vector length.
c
c***references
c***routines called    (none)
c***end prologue       gvadd
c
      real*8 v(nv,n),w(nw,n),x(nx,n)
c
      do 1 i=1,n
         v(1,i)=w(1,i)+x(1,i)
    1 continue
c
      return
      end
