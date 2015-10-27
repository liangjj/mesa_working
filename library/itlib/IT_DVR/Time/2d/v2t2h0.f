*deck v2t2h0.f
c***begin prologue     v2t2h0
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
c***end prologue       v2t2h0
      subroutine v2t2h0(vecin,vecout,tmp,u1,u2,n1,n2,nt,nc,nvc)
      implicit integer (a-z)
      real*8 u1, u2, vecout, vecin, tmp
      dimension u1(n1,n1), u2(n2,n2)
      dimension vecin(n2,n1,nt,nc,2,nvc), vecout(n2,n1,nt,nc,2,nvc)
      dimension tmp(n2,n1,nt,nc,2,nvc)
      common/io/inp, iout
      call ebtc(tmp,u2,vecin,n2,n2,n1*nt*nc*2*nvc)
      do 10 i=1,nt
         do 20 ic=1,nc
            do 30 j=1,nvc
               call ebc(vecout(1,1,i,ic,1,j),
     1                  tmp(1,1,i,ic,1,j),u1,n2,n1,n1) 
               call ebc(vecout(1,1,i,ic,2,j),
     1                  tmp(1,1,i,ic,2,j),u1,n2,n1,n1) 
 30         continue   
 20      continue
 10   continue    
      return
      end       
