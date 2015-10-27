*deck f1s.f 
c***begin prologue     f1s
c***date written       020303   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   
c***author             schneider, b. i.(nsf)
c***source             hcoll
c***purpose            normalized 1s orbital 
c***                   
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       f1s
      subroutine f1s(psi_1s,vpsi_1s,r,n)
c
      implicit integer (a-z)
      real*8 psi_1s, vpsi_1s, r
      dimension psi_1s(n), r(n), vpsi_1s(n)
      common/io/inp, iout    
      do 10 i=1,n
         psi_1s(i)=2.d0*r(i)*exp(-r(i))
         vpsi_1s(i) = psi_1s(i)/r(i)
 10   continue   
      return
      end


















