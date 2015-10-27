*deck ebcxx
      subroutine ebcxx(a,b,c,ni,nk,nj,na,nb,nc)
c***begin prologue     ebcxx
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply,
c***author             saxe, paul (lanl)
c***source             @(#)ebcxx.f        5.1   7/29/88
c***purpose            vectorized matrix operation:  a=b*c .
c***description
c                      call ebcxx(a,b,c,ni,nk,nj,na,nb,nc)
c                        a       output matrix, (na,nj).
c                        b       input matrix, (nb,nk).
c                        c       input matrix, (nc,nj).
c
c***references
c***routines called    sgemm(clams)
c***end prologue       ebcx
      implicit real*8 (a-h,o-z)
c
      parameter (zero=0.0d+0,one=1.0d+0)
      dimension a(na,nj),b(nb,nk),c(nc,nj)
c
      call sgemm('n','n',ni,nj,nk,one,b,nb,c,nc,zero,a,na)
c      call sgmm(ni,nk,nj,b,nb,c,nc,a,na,0,1)
c      call rzero(a,ni*nj)
c      call mxmb(b,1,ni,c,1,nk,a,1,ni,ni,nk,nj)
c
      return
      end
