*deck ebtcxx
      subroutine ebtcxx(a,b,c,ni,nk,nj,na,nb,nc)
c***begin prologue     ebtcxx
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply, subtract
c***author             saxe, paul (lanl)
c***source             @(#)ebtcxx.f      5.1   7/29/88
c***purpose                                             t
c                      vectorized matrix operation:  a=b *c .
c***description
c                      call ebtcxx(a,b,c,ni,nk,nj,na,nb,nc)
c                        b         input matrix, (nk,ni).
c                        c         input matrix, (nk,nj).
c
c***references
c***routines called    sgemm(clams)
c***end prologue       ebtcxx
      implicit real*8 (a-h,o-z)
c
      parameter (zero=0.0d+0,one=1.0d+0)
      dimension a(na,nj),b(nb,ni),c(nc,nj)
c
c
      call sgemm('t','n',ni,nj,nk,one,b,nb,c,nc,zero,a,na)
c      call sgmm(ni,nk,nj,b,nb,c,nc,a,na,4,1)
c      call rzero(a,ni*nj)
c      call mxmb(b,ni,1,c,1,nk,a,1,ni,ni,nk,nj)
c
      return
      end
