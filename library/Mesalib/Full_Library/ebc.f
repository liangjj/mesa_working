*deck @(#)ebc.f	5.1  11/6/94
      subroutine ebc(a,b,c,ni,nk,nj)
c***begin prologue     ebc
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply,
c***author             saxe, paul (lanl)
c***source             @(#)ebc.f	5.1   11/6/94
c***purpose            vectorized matrix operation:  a=b*c .
c***description
c                      call ebc(a,b,c,ni,nk,nj)
c                        a       output matrix, (ni,nj).
c                        b       input matrix, (ni,nk).
c                        c       input matrix, (nk,nj).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       ebc
      implicit integer(a-z)
c
      real*8 a(ni,nj),b(ni,nk),c(nk,nj)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
c     call sgmm(ni,nk,nj,b,ni,c,nk,a,ni,0,1)
c     call rzero(a,ni*nj)
c     call mxmb(b,1,ni,c,1,nk,a,1,ni,ni,nk,nj)
c
c     call mxma(b,1,ni,c,1,nk,a,1,ni,ni,nk,nj)
      call sgemm('n','n',ni,nj,nk,one,b,ni,c,nk,zero,a,ni)
c
      return
      end
