*deck v3t2h0.f
c***begin prologue     v3t2h0
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
c***end prologue       v3t2h0
      subroutine v3t2h0(vecin,vecout,tmp,u1,u2,u3,n1,n2,n3,nt,nc,nvc)
      implicit integer (a-z)
      real*8 u1, u2, u3, vecout, vecin, tmp
      dimension u1(n1,n1), u2(n2,n2), u3(n3,n3)
      dimension vecin(n3,n2,n1,nt,nc,2,nvc)
      dimension vecout(n3,n2,n1,nt,nc,2,nvc)
      dimension tmp(n3,n2,n1,nt,nc,2,nvc)
      common/io/inp, iout
      call ebtc(vecout,u3,vecin,n3,n3,n2*n1*nt*nc*2*nvc)
      do 10 i=1,n1
         do 20 j=1,nt
            do 30 ic=1,nc
               do 40 k=1,nvc
                  call ebc(tmp(1,1,i,j,ic,1,k),
     1                  vecout(1,1,i,j,ic,1,k),u2,n3,n2,n2) 
                  call ebc(tmp(1,1,i,j,ic,2,k),
     1                  vecout(1,1,i,j,ic,2,k),u2,n3,n2,n2) 
 40            continue   
 30         continue
 20      continue   
 10   continue    
      do 50 i=1,nt
         do 60 ic=1,nc
            do 70 j=1,nvc
               call ebc(vecout(1,1,1,i,ic,1,j),
     1                     tmp(1,1,1,i,ic,1,j),u1,n3*n2,n1,n1)
               call ebc(vecout(1,1,1,i,ic,2,j),
     1                     tmp(1,1,1,i,ic,2,j),u1,n3*n2,n1,n1)
 70         continue   
 60      continue
 50   continue   
      return
      end       
