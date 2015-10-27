*deck @(#)atpbct.f	5.1  11/6/94
      subroutine atpbct(a,b,c,ni,nk,nj)
c***begin prologue     atpbct
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply, add
c***author             saxe, paul (lanl)
c***source             @(#)atpbct.f	5.1   11/6/94
c***purpose                                             t    t
c                      vectorized matrix operation:  a=a +b*c .
c***description
c                      call atpbct(a,b,c,ni,nk,nj)
c                        a       output matrix, (nj,ni).
c                        b       input matrix, (ni,nk).
c                        c       input matrix, (nj,nk).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       atpbct
c
      implicit real*8 (a-h,o-z)
c
      real*8 a(nj,ni),b(ni,nk),c(nj,nk)
c
      call sgmm(ni,nk,nj,b,ni,c,nj,a,nj,3,2)
c
      return
      end
