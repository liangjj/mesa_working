*deck ambtcx
      subroutine ambtcx(a,na,b,nb,c,nc,ni,nk,nj)
c***begin prologue     ambtcx
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           matrix, multiply, subtract
c***author             schneider, barry (lanl)
c***source             @(#)ambtcx.f     5.1   8/04/89
c***purpose                                               t
c                      vectorized matrix operation:  a=a-b *c .
c***description
c                      call ambtcx(a,na,b,nb,c,nc,ni,nk,nj)
c                        a       output matrix, (ni,nj).
c                        b       input matrix, (nk,ni).
c                        c       input matrix, (nk,nj).
c
c***references
c***routines called    sgemm(clams)
c***end prologue       ambtcx
      implicit real*8 (a-h,o-z)
c
      dimension a(na,nj),b(nb,ni),c(nc,nj)
c
      parameter (zero=0.0d+0,one=1.0d+0)
c
      call sgemm('t','n',ni,nj,nk,-one,b,nb,c,nc,one,a,na)
c      call sgmm(ni,nk,nj,b,nb,c,nc,a,na,4,-2)
c
      return
      end
