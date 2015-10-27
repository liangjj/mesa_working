*deck cehbtc
c***begin prologue     cehbtc
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           complex matrix multiply
c***author             schneider, barry (lanl)
c***source             mylib
c***                                         +
c***purpose            matrix multiply  a = b * c
c***description        vectorized matrix multiply
c***                   for complex a, b and c
c***                  
c***references         none
c
c***routines called    czero(mylib)
c***end prologue       cehbtc
      subroutine cehbtc(a,b,c,ni,nk,nj)
      implicit integer (a-z)
      complex *16 a, b, c, zero, one
      dimension a(ni,nj), b(nk,ni), c(nk,nj)
      zero=cmplx(0.d0,0.d0)
      one=cmplx(1.d0,0.d0)
      call cgemm('c','n',ni,nj,nk,one,b,nk,c,nk,zero,a,ni)
      return
      end
