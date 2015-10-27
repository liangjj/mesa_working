*deck mkrhs.f
c***begin prologue     mkrhs
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            initial state wavefunction for inhomogeneous, 
c***                   time-dependent hamiltonian.
c***                   
c***references         
c
c***routines called    
c***end prologue       mkrhs
      subroutine mkrhs(rhs,n,m)
      implicit integer (a-z)
      real*8 rhs
      dimension rhs(n,m) 
      common/io/inp, iout
      call rzero(rhs,n*m)
      do 10 i=1,m
         rhs(i,i)=1.d0
   10 continue      
      return
      end       
