*deck mkrhs.f
c***begin prologue     mkrhs
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            
c***                   
c***references         
c
c***routines called    
c***end prologue       mkrhs
      subroutine mkrhs(rhs,n,nrhs)
      implicit integer (a-z)
      real*8 rhs
      dimension rhs(n,nrhs) 
      common/io/inp, iout
      call rzero(rhs,n*nrhs)
      do 10 i=1,nrhs
         rhs(i,i)=1.d0
 10   continue   
      return
      end       
