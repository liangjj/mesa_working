*deck rfive
      subroutine rfive(wt,stp)
c***begin prologue     rfive
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            newton cotes five point quadrature 
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       rfive
c
      implicit integer (a-z)
      real*8 wt, stp, d
      dimension wt(5,4)
      data d / 720.d0 /
      common /io/ inp, iout
      wt(1,1)=251.d0/d
      wt(2,1)=646.d0/d
      wt(3,1)=-264.d0/d
      wt(4,1)=106.d0/d
      wt(5,1)=-19.d0/d
      wt(1,2)=-19.d0/d
      wt(2,2)=346.d0/d
      wt(3,2)=456.d0/d
      wt(4,2)=-74.d0/d
      wt(5,2)=11.d0/d
      wt(1,3)=11.d0/d
      wt(2,3)=-74.d0/d
      wt(3,3)=456.d0/d
      wt(4,3)=346.d0/d
      wt(5,3)=-19.d0/d
      wt(1,4)=-19.d0/d
      wt(2,4)=106.d0/d
      wt(3,4)=-264.d0/d
      wt(4,4)=646.d0/d
      wt(5,4)=251.d0/d
      do 10 i=1,5
         do 20 j=1,4
            wt(i,j)=wt(i,j)*stp
   20    continue
   10 continue          
      return
      end















