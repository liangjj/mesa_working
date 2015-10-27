*deck int_9
      subroutine int_9(wt,stp)
c***begin prologue     int_9
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            newton cotes nine point quadrature 
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       int_9
c
      implicit integer (a-z)
      real*8 wt, stp, fac
      dimension wt(9)
      common /io/ inp, iout
      fac=4.d0*stp/14175.d0
      wt(1)=fac*989.d0
      wt(2)=fac*5888.d0
      wt(3)=-fac*928.d0
      wt(4)=fac*10496.d0
      wt(5)=-fac*454.d0
      wt(6)=wt(4)
      wt(7)=wt(3)
      wt(8)=wt(2)
      wt(9)=wt(1)
      return
      end















