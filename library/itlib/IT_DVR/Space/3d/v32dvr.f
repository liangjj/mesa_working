*deck v32dvr.f
c***begin prologue     v32dvr
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
c***end prologue       v32dvr
      subroutine v32dvr(vecin,vecout,tmp,u1,u2,u3,n1,n2,n3,nc,nvc)
      implicit integer (a-z)
      real*8 vecin, vecout, tmp, u1, u2, u3
      dimension vecin(n3,n2,n1,nc,nvc), vecout(n3,n2,n1,nc,nvc)
      dimension tmp(n3,n2,n1,nc,nvc)
      dimension u1(n1,n1), u2(n2,n2), u3(n3,n3)
      common/io/inp, iout
      call ebc(vecout,u3,vecin,n3,n3,n2*n1*nc*nvc)
      do 10 i=1,n1
         do 20 ic=1,nc
            do 30 j=1,nvc
               call ebct(tmp(1,1,i,ic,j),
     1                vecout(1,1,i,ic,j),u2,n3,n2,n2)  
 30         continue
 20      continue   
 10   continue    
      do 40 ic=1,nc
         do 50 j=1,nvc
            call ebct(vecout(1,1,1,ic,j),
     1                   tmp(1,1,1,ic,j),u1,n3*n2,n1,n1)
 50      continue
 40   continue   
      return
      end       
