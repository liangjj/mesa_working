*deck grid.f 
c***begin prologue     grid
c***date written       020303   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   
c***author             schneider, b. i.(nsf)
c***source             hcoll
c***purpose            driver electron + h atom collisions
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       grid
      subroutine grid(r,stp,n)
c
      implicit integer (a-z)
      real*8 r, stp
      dimension r(n)
      common/io/inp, iout      
      r(1)=1.d-10
      do 10 i=2,n
         r(i) =r(i-1) + stp
 10   continue   
      return
      end


















