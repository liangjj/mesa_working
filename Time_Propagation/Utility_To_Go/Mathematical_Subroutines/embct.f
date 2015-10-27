*deck @(#)embct.f	5.1  11/6/94
      subroutine embct(a,b,c,ni,nk,nj)
c***begin prologue     embct
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply,
c***author             saxe, paul (lanl)
c***source             @(#)embct.f	5.1   11/6/94
c***purpose                                                t
c                      vectorized matrix operation:  a=-b*c .
c***description
c                      call embct(a,b,c,ni,nk,nj)
c                        a       output matrix, (ni,nj).
c                        b       input matrix, (ni,nk).
c                        c       input matrix, (nj,nk).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       embct
      implicit  integer(a-z)
c
      real*8 a(ni,nj),b(ni,nk),c(nj,nk)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
c     call sgmm(ni,nk,nj,b,ni,c,nj,a,ni,2,-1)
      call sgemm('n','t',ni,nj,nk,-one,b,ni,c,nj,zero,a,ni)
c
      return
      end
