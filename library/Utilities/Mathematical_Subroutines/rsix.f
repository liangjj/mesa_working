*deck rsix
      subroutine rsix(wt,stp)
c***begin prologue     rsix
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            newton cotes six point quadrature 
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       rsix
c
      implicit integer (a-z)
      real*8 wt, stp, d
      dimension wt(6,5)
      data d / 1440.d0 /
      common /io/ inp, iout
      wt(1,1)=475.d0/d
      wt(2,1)=1427.d0/d
      wt(3,1)=-798.d0/d
      wt(4,1)=482.d0/d
      wt(5,1)=-173.d0/d
      wt(6,1)=27.d0/d
      wt(1,2)=-27.d0/d
      wt(2,2)=637.d0/d
      wt(3,2)=1022.d0/d
      wt(4,2)=-258.d0/d
      wt(5,2)=77.d0/d
      wt(6,2)=-11.d0/d
      wt(1,3)=11.d0/d
      wt(2,3)=-93.d0/d
      wt(3,3)=802.d0/d
      wt(4,3)=802.d0/d
      wt(5,3)=-93.d0/d
      wt(6,3)=11.d0/d
      wt(1,4)=-11.d0/d
      wt(2,4)=77.d0/d
      wt(3,4)=-258.d0/d
      wt(4,4)=1022.d0/d
      wt(5,4)=637.d0/d
      wt(6,4)=-27.d0/d
      wt(1,5)=27.d0/d
      wt(2,5)=-173.d0/d
      wt(3,5)=482.d0/d
      wt(4,5)=-798.d0/d
      wt(5,5)=1427.d0/d
      wt(6,5)=475.d0/d
      do 10 i=1,6
         do 20 j=1,5
            wt(i,j)=wt(i,j)*stp
   20    continue
   10 continue          
      return
      end















