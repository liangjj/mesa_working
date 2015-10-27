*deck @(#)ambc.f	5.1  11/6/94
      subroutine ambc(a,b,c,ni,nk,nj)
c***begin prologue     ambc
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply, subtract
c***author             saxe, paul (lanl)
c***source             @(#)ambc.f	5.1   11/6/94
c***purpose            vectorized matrix operation: a=a-b*c .
c***description
c                      call ambc(a,b,c,ni,nk,nj)
c                        a       output matrix, (ni,nj).
c                        b       input matrix, (ni,nk).
c                        c       input matrix, (nk,nj).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       ambc
      implicit integer(a-z)
c
      real*8 a(ni,nj),b(ni,nk),c(nk,nj)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
      call sgemm('n','n',ni,nj,nk,-one,b,ni,c,nk,one,a,ni)
c     call sgmm(ni,nk,nj,b,ni,c,nk,a,ni,0,-2)
c     call mxmbn(b,1,ni,c,1,nk,a,1,ni,ni,nk,nj)
c
      return
      end
