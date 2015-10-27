*deck int_7
      subroutine int_7(wt,stp)
c***begin prologue     int_7
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            newton cotes seven point quadrature 
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       int_7
c
      implicit integer (a-z)
      real*8 wt, stp, fac
      dimension wt(7)
      common /io/ inp, iout
      fac=stp/144.d0
      wt(1)=fac*41.d0
      wt(2)=fac*216.d0
      wt(3)=fac*27.d0
      wt(4)=fac*272.d0
      wt(5)=wt(3)
      wt(6)=wt(2)
      wt(7)=wt(1)
      return
      end















