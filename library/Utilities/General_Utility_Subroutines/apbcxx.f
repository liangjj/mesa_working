*deck apbcxx
      subroutine apbcxx(a,b,c,ni,nk,nj,na,nb,nc)
c***begin prologue     apbcxx
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply, subtract
c***author             saxe, paul (lanl)
c***source             @(#)apbcx.f       5.1   7/29/88
c***purpose            vectorized matrix operation: a=a+b*c .
c***description
c                      call apbcxx(a,b,c,ni,nk,nj,na,nb,nc)
c                        a       output matrix, (na,nj).
c                        b       input matrix, (nb,nk).
c                        c       input matrix, (nc,nj).
c
c***references
c***routines called    sgemm(clams)
c***end prologue       apbcxx
      implicit real*8 (a-h,o-z)
      parameter ( zero=0.d0,one=1.d0)
c
      dimension a(na,nj),b(nb,nk),c(nc,nj)
c
      call sgemm('n','n',ni,nj,nk,one,b,nb,c,nc,one,a,na)
c     call sgmm(ni,nk,nj,b,nb,c,nc,a,na,0,-2)
c      call mxmbn(b,1,ni,c,1,nk,a,1,ni,ni,nk,nj)
c
      return
      end
