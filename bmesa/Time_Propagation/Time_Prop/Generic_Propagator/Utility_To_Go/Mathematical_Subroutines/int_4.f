*deck int_4
      subroutine int_4(wt,stp)
c***begin prologue     int_4
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
c***end prologue       int_4
c
      implicit integer (a-z)
      real*8 wt, stp
      dimension wt(4)
      common /io/ inp, iout
      wt(1)=3.d0*stp/8.d0
      wt(2)=3.d0*wt(1)
      wt(3)=wt(2)
      wt(4)=wt(1)
      return
      end















