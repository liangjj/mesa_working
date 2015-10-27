*deck @(#)ebtc.f	5.1  11/6/94
      subroutine ebtc(a,b,c,ni,nk,nj)
c***begin prologue     ebtc
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply, subtract
c***author             saxe, paul (lanl)
c***source             @(#)ebtc.f	5.1   11/6/94
c***purpose                                             t
c                      vectorized matrix operation:  a=b *c .
c***description
c                      call ebtc(a,b,c,ni,nk,nj)
c                        b       input matrix, (nk,ni).
c                        c       input matrix, (nk,nj).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       ebtc
      implicit integer(a-z)
c
      real*8 a(ni,nj),b(nk,ni),c(nk,nj)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
c     call sgmm(ni,nk,nj,b,nk,c,nk,a,ni,4,1)
c     call rzero(a,ni*nj)
c     call mxmb(b,nk,1,c,1,nk,a,1,ni,ni,nk,nj)
c
c     call mxma(b,nk,1,c,1,nk,a,1,ni,ni,nk,nj)
      call sgemm('t','n',ni,nj,nk,one,b,nk,c,nk,zero,a,ni)
c
      return
      end
