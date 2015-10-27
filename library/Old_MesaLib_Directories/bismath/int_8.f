*deck int_8
      subroutine int_8(wt,stp)
c***begin prologue     int_8
c***date written       920324   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            newton cotes eight point quadrature 
c***
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       int_8
c
      implicit integer (a-z)
      real*8 wt, stp, fac
      dimension wt(8)
      common /io/ inp, iout
      fac=7.d0*stp/17280.d0
      wt(1)=fac*751.d0
      wt(2)=fac*3577.d0
      wt(3)=fac*1323.d0
      wt(4)=fac*2989.d0
      wt(5)=wt(4)
      wt(6)=wt(3)
      wt(7)=wt(2)
      wt(8)=wt(1)
      return
      end















