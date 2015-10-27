*deck h32dvr.f
c***begin prologue     h32dvr
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            transform a 3-D vector from the dvr to  
c***                   the H0 representation.
c***                                      
c***references         
c
c***routines called    
c***end prologue       h32dvr
      subroutine h32dvr(vecin,vecout,tmp,u1,u2,u3,n1,n2,n3,nvc)
      implicit integer (a-z)
      real*8 vecin, vecout, tmp, u1, u2, u3
      dimension vecin(n3,n2,n1,nvc), vecout(n3,n2,n1,nvc)
      dimension tmp(n3,n2,n1,nvc)
      dimension u1(n1,n1), u2(n2,n2), u3(n3,n3)
      common/io/inp, iout
      call ebc(vecout,u3,vecin,n3,n3,n2*n1*nvc)
      do 10 i=1,n1
         do 20 j=1,nvc
            call ebct(tmp(1,1,i,j),vecout(1,1,i,j),u2,n3,n2,n2)  
 20      continue   
 10   continue    
      do 30 i=1,nvc
         call ebct(vecout(1,1,1,i),tmp(1,1,1,i),u1,n3*n2,n1,n1)
 30   continue   
      return
      end       
