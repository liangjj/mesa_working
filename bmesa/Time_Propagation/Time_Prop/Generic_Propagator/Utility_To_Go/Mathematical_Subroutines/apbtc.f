*deck @(#)apbtc.f	5.1  11/6/94
      subroutine apbtc(a,b,c,ni,nk,nj)
c***begin prologue     apbtc
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply, addt
c***author             saxe, paul (lanl)
c***source             @(#)apbtc.f	5.1   11/6/94
c***purpose                                               t
c                      vectorized matrix operation:  a=a+b *c .
c***description
c                      call apbtc(a,b,c,ni,nk,nj)
c                        a       output matrix, (ni,nj).
c                        b       input matrix, (nk,ni).
c                        c       input matrix, (nk,nj).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       apbtc
      implicit  integer(a-z)
c
      real*8 a(ni,nj),b(nk,ni),c(nk,nj)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
      call sgemm('t','n',ni,nj,nk,one,b,nk,c,nk,one,a,ni)
c     call sgmm(ni,nk,nj,b,nk,c,nk,a,ni,4,2)
c
      return
      end
