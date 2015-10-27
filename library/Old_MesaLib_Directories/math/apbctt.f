*deck @(#)apbctt.f	5.1  11/6/94
      subroutine apbctt(a,b,c,ni,nk,nj)
c***begin prologue     apbctt
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply, add
c***author             saxe, paul (lanl)
c***source             @(#)apbctt.f	5.1   11/6/94
c***purpose                                             t  t  t
c                      vectorized matrix operation:  a=a +b *c .
c***description
c                      call apbctt(a,b,c,ni,nk,nj)
c                        a       output matrix, (nj,ni).
c                        b       input matrix, (nk,ni).
c                        c       input matrix, (nj,nk).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       apbctt
      implicit integer(a-z)
c
      real*8 a(nj,ni),b(nk,ni),c(nj,nk)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
      call sgmm(ni,nk,nj,b,nk,c,nj,a,nj,7,2)
c
      return
      end
