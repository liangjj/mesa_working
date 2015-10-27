*deck cambct
c***begin prologue     cambct
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           complex matrix multiply
c***author             schneider, barry (lanl)
c***source             mylib
c***                                                 t
c***purpose            matrix multiply  a = a - b * c
c***description        vectorized matrix multiply
c***                   for complex a, b and c
c***                  
c***references         none
c
c***routines called    czero(mylib)
c***end prologue       cambct 
      subroutine cambct(a,b,c,ni,nk,nj)
      implicit integer (a-z)
      complex *16 a,b,c,zero, one
      dimension a(ni,nj), b(ni,nk), c(nj,nk)
      zero=cmplx(0.d0,0.d0)
      one=cmplx(1.d0,0.d0)
      call cgemm('n','t',ni,nj,nk,-one,b,ni,c,nj,one,a,ni)     
      return
      end
