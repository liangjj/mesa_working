*deck int_6
      subroutine int_6(wt,stp)
c***begin prologue     int_6
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
c***end prologue       int_6
c
      implicit integer (a-z)
      real*8 wt, stp, fac
      dimension wt(6)
      common /io/ inp, iout
      fac=5.d0*stp/288.d0
      wt(1)=fac*19.d0
      wt(2)=fac*75.d0
      wt(3)=fac*50.d0
      wt(4)=wt(3)
      wt(5)=wt(2)
      wt(6)=wt(1)
      return
      end















