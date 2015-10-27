*deck guess.f
c***begin prologue     guess
c***date written       960615   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hamiltonian
c***author             schneider, barry (nsf)
c***source             3-dim
c***purpose            guess vectors based on separable hamiltonian
c***                   for davidson routine.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       guess
      subroutine guess(vec,eig,n,nroots)
      implicit integer (a-z)
      real*8 vec, eig
      dimension vec(n,nroots), eig(nroots)
      common/io/inp, iout
      call rzero(vec,n*nroots)
      do 10 i=1,nroots
         eig(i)=1.d0 
         vec(i,i)=1.0d0
 10   continue   
      return
      end       






