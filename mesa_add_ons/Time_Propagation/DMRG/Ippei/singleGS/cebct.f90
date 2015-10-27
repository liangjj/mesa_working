!deck cebct
!***begin prologue     cebct
!***date written       880423   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           complex matrix multiply
!***author             schneider, barry (lanl)
!***source             mylib
!                                                t
!***purpose            matrix multiply  a = b * c
!***description        vectorized matrix multiply
!***                   for complex a, b and c
!***                  
!***references         none
!
!***routines called    czero(mylib)
!***end prologue       cebct
      subroutine cebct(a,b,c,ni,nk,nj)
      implicit none
      integer                           :: ni, nj, nk
      complex*16, dimension(ni,nj)      :: a
      complex*16, dimension(ni,nk)      :: b
      complex*16, dimension(nj,nk)      :: c
      complex*16                        :: zero=(0.d0,0.d0), one=(1.d0,0.d0)
       call zgemm('n','t',ni,nj,nk,one,b,ni,c,nj,zero,a,ni)
      return
      end
