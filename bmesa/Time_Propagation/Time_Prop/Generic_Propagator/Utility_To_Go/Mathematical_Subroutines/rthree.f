*deck rthree
      subroutine rthree(wt,stp)
c***begin prologue     rthree
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            newton cotes for three point quadrature
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       rthree
c
      implicit integer (a-z)
      real*8 wt, stp, d
      dimension wt(3,2)
      data d / 12.d0 /
      common /io/ inp, iout
      wt(1,1)=5.d0/d
      wt(2,1)=8.d0/d
      wt(3,1)=-1.d0/d
      wt(1,2)=wt(3,1)
      wt(2,2)=wt(2,1)
      wt(3,2)=wt(1,1)
      do 10 i=1,3
         do 20 j=1,2
            wt(i,j)=wt(i,j)*stp
   20    continue
   10 continue          
      return
      end















