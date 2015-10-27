*deck rtwo
      subroutine rtwo(wt,stp)
c***begin prologue     rtwo
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            newton cotes for two point quadrature 
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       rtwo
c
      implicit integer (a-z)
      real*8 wt, stp
      dimension wt(2,1)
      common /io/ inp, iout
      wt(1,1)=stp*.5d0
      wt(2,1)=wt(1,1) 
      return
      end















