*deck local.f 
c***begin prologue     local
c***date written       020303   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   
c***author             schneider, b. i.(nsf)
c***source             hcoll
c***purpose            normalized 1s orbital, 
c***                   local potential for h + e collisions 
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       local
      subroutine local(v_loc,r,n)
c
      implicit integer (a-z)
      real*8 v_loc, r
      dimension v_loc(n), r(n)
      common/io/inp, iout    
      do 10 i=1,n
c         v_loc(i) =  ( 1.d0 + 1.d0/r(i) )*exp(-2.d0*r(i))
          v_loc(i) = exp(-r(i))
 10   continue   
      return
      end


















