*deck mkpsi0.f
c***begin prologue     mkpsi0
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
c***end prologue       mkpsi0
      subroutine mkpsi0(v,rhs,n,state)
      implicit integer (a-z)
      real*8 v, rhs
      dimension v(n,n), rhs(n) 
      common/io/inp, iout
      call iosys('read real "H0 eigenfunctions" from lamdat',
     1            n*n,v,0,' ')
      call copy(v(1,state+1),rhs,n)
      return
      end       
