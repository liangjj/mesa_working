c $Header: frms.f,v 1.1 92/12/31 14:43:12 bis Exp $
*deck frms.f
c***begin prologue     frms
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           rms
c***author             schneider, barry (nsf)
c***source             m6201
c***purpose            rms error of function 
c***references         none
c
c***routines called
c***end prologue       frms
      subroutine frms (fa,fb,rms,n)
      implicit integer (a-z)
      real*8 fa, fb, rms
      dimension fa(n), fb(n)
      rms=0.d0
      do 10 i=1,n
         rms=rms+(fa(i)-fb(i))*(fa(i)-fb(i))
   10 continue
      rms=sqrt(rms)/n
      return
      end
