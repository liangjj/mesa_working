*deck @(#)apbct.f	5.1  11/6/94
      subroutine apbct(a,b,c,ni,nk,nj)
c***begin prologue     apbct
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply, add
c***author             saxe, paul (lanl)
c***source             @(#)apbct.f	5.1   11/6/94
c***purpose                                                t
c                      vectorized matrix operation: a=a+b*c .
c***description
c                      call apbct(a,b,c,ni,nk,nj)
c                        a       output matrix, (ni,nj).
c                        b       input matrix, (ni,nk).
c                        c       input matrix, (nj,nk).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       apbct
      implicit  integer(a-z)
c
      real*8 a(ni,nj),b(ni,nk),c(nj,nk)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
      call sgemm('n','t',ni,nj,nk,one,b,ni,c,nj,one,a,ni)
c     call sgmm(ni,nk,nj,b,ni,c,nj,a,ni,2,2)
c
      return
      end
