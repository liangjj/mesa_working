*deck @(#)ebtct.f	5.1  11/6/94
      subroutine ebtct(a,b,c,ni,nk,nj)
c***begin prologue     ebtct
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply,
c***author             saxe, paul (lanl)
c***source             @(#)ebtct.f	5.1   11/6/94
c***purpose                                             t  t
c                      vectorized matrix operation:  a=b *c .
c***description
c                      call ebtct(a,b,c,ni,nk,nj)
c                        a       output matrix, (ni,nj).
c                        b       input matrix, (nk,ni).
c                        c       input matrix, (nj,nk).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       ebtct
      implicit  integer(a-z)
c
      real*8 a(ni,nj),b(nk,ni),c(nj,nk)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
c     call sgmm(ni,nk,nj,b,nk,c,nj,a,ni,6,1)
c     call rzero(a,ni*nj)
c     call mxmb(b,nk,1,c,nj,1,a,1,ni,ni,nk,nj)
c
c     call mxma(b,nk,1,c,nj,1,a,1,ni,ni,nk,nj)
      call sgemm('t','t',ni,nj,nk,one,b,nk,c,nj,zero,a,ni)
c
      return
      end
