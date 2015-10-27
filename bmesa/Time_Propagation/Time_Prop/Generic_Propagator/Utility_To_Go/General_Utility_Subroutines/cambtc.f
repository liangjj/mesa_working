*deck cambtc
c***begin prologue     cambtc
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           complex matrix multiply
c***author             schneider, barry (lanl)
c***source             mylib
c***                                             t
c***purpose            matrix multiply  a = a - b * c
c***description        vectorized matrix multiply
c***                   for complex a, b and c
c***                  
c***references         none
c
c***routines called    czero(mylib)
c***end prologue       cambtc
      subroutine cambtc(a,b,c,ni,nk,nj)
      implicit integer (a-z)
      complex *16 a, b, c, zero, one
      dimension a(ni,nj), b(nk,ni), c(nk,nj)
c
      zero=cmplx(0.d0,0.d0)
      one=cmplx(1.d0,0.d0) 
      call cgemm('t','n',ni,nj,nk,-one,b,nk,c,nk,one,a,ni)
      return
      end
