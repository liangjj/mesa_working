*deck h12h0.f
c***begin prologue     h12h0
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            transform a 1-D vector from the dvr to  
c***                   the H0 representation.
c***                                      
c***references         
c
c***routines called    
c***end prologue       h12h0
      subroutine h12h0(vecin,vecout,u1,n1,nvc)
      implicit integer (a-z)
      real*8 vecout, vecin, u1
      dimension vecin(n1,nvc), vecout(n1,nvc), u1(n1,n1)
      common/io/inp, iout
      call ebtc(vecout,u1,vecin,n1,n1,nvc)
      return
      end       
