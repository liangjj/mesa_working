*deck integ.f 
c***begin prologue     local
c***date written       020303   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   
c***author             schneider, b. i.(nsf)
c***source             hcoll
c***purpose            integration
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       integ
      function integ(f1,f2,h,n)
c
      implicit integer (a-z)
      real*8 integ, f1, f2, h, sdot
      dimension f1(n), f2(n)
      common/io/inp, iout    
      integ= ( f1(1)*f2(1) + f1(n)*f2(n) ) * .5d0
      integ = integ + sdot(n-2,f1(2),1,f2(2),1)
      integ = h*integ
      return
      end


















