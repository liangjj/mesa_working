*deck @(#)atebc.f	5.1  11/6/94
      subroutine atebc(a,b,c,ni,nk,nj)
c***begin prologue     atebc
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply,
c***author             saxe, paul (lanl)
c***source             @(#)atebc.f	5.1   11/6/94
c***purpose                                           t
c                      vectorized matrix operation:  a =b*c.
c***description
c                      call atebc(a,b,c,ni,nk,nj)
c                        a       output matrix,(nj,ni).
c                        b       input matrix, (ni,nk).
c                        c       input matrix, (nk,nj).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       atebc
      implicit integer(a-z)
c
      real*8 a(nj,ni),b(ni,nk),c(nk,nj)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
      call sgmm(ni,nk,nj,b,ni,c,nk,a,nj,1,1)
c
      return
      end
