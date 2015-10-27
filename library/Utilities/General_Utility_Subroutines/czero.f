*deck czero
c***begin prologue     czero
c***date written       880423   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           complex array
c***author             schneider, barry (lanl)
c***source             mylib
c***purpose            zero complex array
c***description        set complex array to zero
c***                  
c***references         none
c
c***routines called    none
c***end prologue       czero
      subroutine czero(array,n)
      dimension array(n)
      complex *16 array
      do 100 i=1,n
  100 array(i) = cmplx(0.d+00,0.d+00)
      return
      end
