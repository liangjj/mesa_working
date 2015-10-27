*deck h22h0.f
c***begin prologue     h22h0
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            transform a 2-D vector from the dvr to  
c***                   the H0 representation.
c***                                      
c***references         
c
c***routines called    
c***end prologue       h22h0
      subroutine h22h0(vecin,vecout,tmp,u1,u2,n1,n2,n,nvc)
      implicit integer (a-z)
      real*8 vecin, vecout, tmp, u1, u2
      dimension vecin(n2,n1,nvc), vecout(n2,n1,nvc)
      dimension tmp(n2,n1,nvc)
      dimension u1(n1,n1), u2(n2,n2)
      common/io/inp, iout
      call ebtc(tmp,u2,vecin,n2,n2,n1*nvc)
      do 10 i=1,nvc
         call ebc(vecout(1,1,i),tmp(1,1,i),u1,n2,n1,n1)  
 10   continue    
      return
      end       
