*deck @(#)atembc.f	5.1  11/6/94
      subroutine atembc(a,b,c,ni,nk,nj)
c***begin prologue     atembc
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply
c***author             saxe, paul (lanl)
c***source             @(#)atembc.f	5.1   11/6/94
c***purpose                                           t
c                      vectorized matrix operation:  a =-b*c.
c***description
c                      call atembc(a,b,c,ni,nk,nj)
c                        b       input matrix, (ni,nk).
c                        c       input matrix, (nk,nj).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       atembc
      implicit real*8 (a-h,o-z)
c
      real*8 a(nj,ni),b(ni,nk),c(nk,nj)
c
      call sgmm(ni,nk,nj,b,ni,c,nk,a,nj,1,-1)
c
      return
      end
