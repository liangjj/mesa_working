*deck @(#)ebct.f	5.1  11/6/94
      subroutine ebct(a,b,c,ni,nk,nj)
c***begin prologue     ebct
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply,
c***author             saxe, paul (lanl)
c***source             @(#)ebct.f	5.1   11/6/94
c***purpose                                               t
c                      vectorized matrix operation:  a=b*c .
c***description
c                      call ebct(a,b,c,ni,nk,nj)
c                        a       output matrix, (ni,nj).
c                        b       input matrix, (ni,nk).
c                        c       input matrix, (nj,nk).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       ebct
      implicit integer(a-z)
c
      real*8 a(ni,nj),b(ni,nk),c(nj,nk)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
c     call sgmm(ni,nk,nj,b,ni,c,nj,a,ni,2,1)
c     call rzero(a,ni*nj)
c     call mxmb(b,1,ni,c,nj,1,a,1,ni,ni,nk,nj)
c
c     call mxma(b,1,ni,c,nj,1,a,1,ni,ni,nk,nj)
      call sgemm('n','t',ni,nj,nk,one,b,ni,c,nj,zero,a,ni)
c
      return
      end
