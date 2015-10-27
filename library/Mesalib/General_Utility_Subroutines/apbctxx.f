*deck @(#)apbctxx.f	5.1  11/6/94
      subroutine apbctxx(a,b,c,ni,nk,nj,na,nb,nc)
c***begin prologue     apbctxx
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply, add
c***author             saxe, paul (lanl)
c***source             @(#)apbct.f	5.1   11/6/94
c***purpose                                                t
c                      vectorized matrix operation: a=a+b*c .
c***description
c                      call apbctxx(a,b,c,ni,nk,nj,na,nb,nc)
c                        a       output matrix, (ni,nj).
c                        b       input matrix, (ni,nk).
c                        c       input matrix, (nj,nk).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       apbctxx
      implicit  integer(a-z)
c
      real*8 a(na,nj), b(nb,nk), c(nc,nk)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
      call sgemm('n','t',ni,nj,nk,one,b,nb,c,nc,one,a,na)
c     call sgmm(ni,nk,nj,b,ni,c,nj,a,ni,2,2)
c
      return
      end
