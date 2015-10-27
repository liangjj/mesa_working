*deck @(#)vscale.f	1.1  11/30/90
      subroutine vscale(v,w,fac,n)
c***begin prologue     vscale
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           vector, scale
c***author             saxe, paul (lanl)
c***source             @(#)vscale.f	1.1   11/30/90
c***purpose                                             
c                      vectorized vector scale:  v=fac*w   .
c***description
c                      call vscale(v,w,fac,n)
c                        w        input vector of length n
c                        v        output vector of length n.
c                        fac      scalar multiplying factor
c                        n        vector length.
c
c***references
c***routines called    (none)
c***end prologue       vscale
      real*8 v(n), w(n), fac
      do 1 i=1,n
         v(i)=fac*w(i)
    1 continue
      return
      end
