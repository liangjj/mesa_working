*deck @(#)ambctxx.f	5.1  11/6/94
      subroutine ambctxx(a,b,c,ni,nk,nj,na,nb,nc)
c***begin prologue     ambctxx
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply, subtract
c***author             saxe, paul (lanl)
c***source             @(#)ambct.f	5.1   11/6/94
c***purpose                                                t
c                      vectorized matrix operation: a=a-b*c .
c***description
c                      call ambctxx(a,b,c,ni,nk,nj,na,nb,nc)
c                        a       output matrix, (ni,nj).
c                        b       input matrix, (ni,nk).
c                        c       input matrix, (nj,nk).
c
c***references
c***routines called    sgmm(clams)
c***end prologue       ambct
      implicit integer(a-z)
c
      real*8 a(na,nj), b(nb,nk), c(nc,nk)
      real*8 zero,one
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
      call sgemm('n','t',ni,nj,nk,-one,b,nb,c,nc,one,a,na)
c     call sgmm(ni,nk,nj,b,nb,c,nc,a,na,2,-2)
c
      return
      end
