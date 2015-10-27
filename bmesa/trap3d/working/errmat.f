*deck errmat.f
c***begin prologue     errmat
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski
c***author             schneider, barry (nsf)
c***source             tstdiis
c***purpose            calculate the error vector from the commutator
c***                   of the fock matrix and density matrix.
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       errmat
      subroutine errmat(err,psi,hpsi,n)
      implicit integer (a-z)
      real*8 err, psi, hpsi, eij
      dimension err(n*n), psi(n), hpsi(n)
      common/io/inp, iout
      ij=0
      do 10 i=1,n
         do 20 j=1,n
            ij = ij + 1
            err(ij) = hpsi(i)*psi(j) - hpsi(j)*psi(i)
 20      continue
 10   continue   
      return 
      end       
