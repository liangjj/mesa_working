*deck rfour
      subroutine rfour(wt,stp)
c***begin prologue     rfour
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            newton cotes for four point quadrature
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       rfour
c
      implicit integer (a-z)
      real*8 wt, stp, d
      dimension wt(4,3)
      data d / 24.d0 /
      common /io/ inp, iout
      wt(1,1)=9.d0/d
      wt(2,1)=19.d0/d
      wt(3,1)=-5.d0/d
      wt(4,1)=1.d0/d
      wt(1,2)=-1.d0/d
      wt(2,2)=13.d0/d
      wt(3,2)=13.d0/d
      wt(4,2)=-1.d0/d
      wt(1,3)=1.d0/d
      wt(2,3)=-5.d0/d
      wt(3,3)=19.d0/d
      wt(4,3)=9.d0/d
      do 10 i=1,4
         do 20 j=1,3
            wt(i,j)=wt(i,j)*stp
   20    continue
   10 continue          
      return
      end















