*deck exact.f
c***begin prologue     exact
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           3-d schroedinger equation
c***author             schneider, barry (nsf)
c***source             trap3d
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       exact
      subroutine exact(rhs,crhs,q,n,mattyp,type)
      implicit integer (a-z)
      real*8 rhs, q
      complex*16 crhs
      character*80 title
      character*(*) type
      character*(*) mattyp
      dimension crhs(n), rhs(n), q(n)
      common/io/inp, iout
      if(mattyp.eq.'complex') then
         if(type.eq.'delta') then
            do 10 i=1,n
               crhs(n) = -2.d0*q(i)
 10         continue               
         elseif(type.eq.'one') then
            do 20 i=1,n
               crhs(i) = q(i)*( q(i) - 2.d0 )
 20         continue
         elseif(type.eq.'x') then
            do 30 i=1,n
               crhs(i) = q(i) * ( q(i)*q(i)/3.d0 - 1.d0 )
 30         continue
         endif                                     
      else
         if(type.eq.'delta') then
            do 40 i=1,n
               rhs(n) = -2.d0*q(i)
 40         continue               
         elseif(type.eq.'one') then
            do 50 i=1,n
               rhs(i) = q(i)*( q(i) - 2.d0 )
 50         continue
         elseif(type.eq.'x') then
            do 60 i=1,n
               rhs(i) = q(i) * ( q(i)*q(i)/3.d0 - 1.d0 )
 60         continue
         endif              
      endif            
      return
      end       
