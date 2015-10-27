!deck gs_energy.f
!***begin prologue     gs_energy
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           eigenvalues
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculates eigenvalues from propagated
!***                   TD wavefunction. 
!***description        Vector_0 is the solution and vector is the
!***                   result of the operation of the Hamiltonian
!***                   on the vector_0.  This operation must
!***                   be performed before entry into gs_energy 
!***references
!***routines called    iosys, util and mdutil
!***end prologue       gs_energy
  SUBROUTINE gs_energy(vector_0,vector,lambda)
  USE arnoldi_global
  USE dvr_shared  
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)                 :: vector_0, vector
  REAL*8                                 :: ddot
  REAL*8                                 :: numerator, denominator
  REAL*8                                 :: lambda
!
  call iosys('read real solution from bec',n3d,vector_0,0,' ')
  numerator   = ddot(n3d,vector_0,1,vector,1)    
  denominator = ddot(n3d,vector,1,vector,1)       
  lambda = numerator / denominator
  lambda = log(lambda)/deltat
END SUBROUTINE gs_energy
!***********************************************************************
