*deck precon.f
c***begin prologue     precon
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative solver
c***author             schneider, barry (nsf)
c***source             
c***purpose            set up hamiltonian and right hand side in a 
c***                   preconditioned form based on the diagonal as a
c***                   preconditioner.
c***references         
c
c***routines called    
c***end prologue       precon
      subroutine precon(ham,diag,rhs,trial,energy,n,nrhs,m)
      implicit integer (a-z)
      real*8 ham, diag, rhs, trial, energy
      dimension ham(n,n), diag(n), rhs(n,nrhs), trial(n,m)
      common/io/inp, iout
      do 10 i=1,n
         diag(i) = ham(i,i) - energy
         ham(i,i) = 0.d0
 10   continue
      call vimmul(diag,ham,ham,n,n)
      call vimmul(diag,trial,trial,n,m)
      call vimmul(diag,rhs,rhs,n,nrhs)
      return
      end       
