*deck @(#)vinv.f	1.2  10/27/94
      subroutine vinv(v,w,n)
c***begin prologue     vinv
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, invert
c***author             saxe, paul (lanl)
c***source             @(#)vinv.f	1.2   10/27/94
c***purpose                                             -1
c                      vectorized vector inversion:  v=w   .
c***description
c                      call vinv(v,w,n)
c                        w        input vector of length n
c                        v        output vector of length n.
c                        n        vector length.
c
c***references
c***routines called    (none)
c***end prologue       vinv
      real*8 v(n),w(n)
      do 1 i=1,n
         v(i)=1.0d+00/w(i)
    1 continue
      return
      end
