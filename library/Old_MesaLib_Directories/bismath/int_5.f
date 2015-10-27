*deck int_5
      subroutine int_5(wt,stp)
c***begin prologue     int_5
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
c***end prologue       int_5
c
      implicit integer (a-z)
      real*8 wt, stp
      dimension wt(5)
      common /io/ inp, iout
      wt(1)=14.d0*stp/45.d0
      wt(2)=64.d0*stp/45.d0
      wt(3)=24.d0*stp/45.d0
      wt(4)=wt(2)
      wt(5)=wt(1)
      return
      end















