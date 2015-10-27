*deck int_3
      subroutine int_3(wt,stp)
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
c***end prologue       int_3
c
      implicit integer (a-z)
      real*8 wt, stp
      dimension wt(3)
      common /io/ inp, iout
      wt(1)=1.d0*stp/3.d0
      wt(2)=4.d0*stp/3.d0
      wt(3)=wt(1)
      return
      end















