*deck @(#)ambtct.f	5.1  11/6/94
      subroutine ambtct(a,b,c,ni,nk,nj)
c***begin prologue     ambtct
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply, subtract
c***author             saxe, paul (lanl)
c***source             @(#)ambtct.f	5.1   11/6/94
c***purpose                                               t  t
c                      vectorized matrix operation:  a=a-b *c .
c***description
c                      call ambtct(a,b,c,ni,nk,nj)
c                        a       output matrix, (ni,nj).
c                        b       input matrix, (nk,ni).
c                        c       input matrix, (nj,nk).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       ambtct
      implicit integer(a-z)
c
      real*8 a(ni,nj),b(nk,ni),c(nj,nk)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
      call sgemm('t','t',ni,nj,nk,-one,b,nk,c,nj,one,a,ni)
c     call sgmm(ni,nk,nj,b,nk,c,nj,a,ni,6,-2)
c
      return
      end
