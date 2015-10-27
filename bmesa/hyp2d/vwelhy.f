*deck vwelhy.f
c***begin prologue     vwelhy
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           potential, hyperspherical
c***author             schneider, barry (nsf)
c***source             
c***purpose            potential well in hyperspherical coordinates
c***                   
c***description        
c***                   
c***references         
c
c***routines called    
c***end prologue       vwelhy
      subroutine vwelhy(v,phi,r,d,len,n)
      implicit integer (a-z)
      real*8 v, r, phi
      real*8 r1, r2, r12, d, len
      dimension d(2,2), len(2,2)
      dimension n(2)
      dimension v(n(1),n(2)), r(n(2)), phi(n(1))
      common/io/inp, iout
      do 10 i=1,n(1)
         ci=cos(phi(i))
         si=sin(phi(i))
         do 20 j=1,n(2)
            r1=r(j)*ci
            r2=r(j)*si
            r12=abs(r2-r1)
            if(r1.le.len(1,1)) then
               v(i,j) = v(i,j) + d(1,1)
            endif      
            if(r2.le.len(2,2)) then
               v(i,j) = v(i,j) + d(2,2)
            endif              
            if(r12.le.len(1,2)) then
               v(i,j) = v(i,j) + d(1,2)
            endif      
 20      continue
 10   continue  
      return
      end       














