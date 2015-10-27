!deck Make_Hamiltonian.f
!***begin prologue     Make_Hamiltonian
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Make_Hamiltonian
!
  SUBROUTINE Make_Hamiltonian(ke,v,eig,eigvec,srf_prm,srf,n)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  INTEGER                                :: i
  INTEGER                                :: info
  REAL*8, DIMENSION(n,n)                 :: ke
  REAL*8, DIMENSION(n)                   :: v
  REAL*8, DIMENSION(n)                   :: eig
  REAL*8, DIMENSION(n,n)                 :: eigvec
  REAL*8                                 :: srf_prm
  REAL*8, DIMENSION(n)                   :: srf
  CHARACTER(LEN=80)                      :: title
  REAL*8, DIMENSION(:), ALLOCATABLE      :: work
  ALLOCATE(work(5*n)) 
  eigvec = ke
  DO i=1,n
     eigvec(i,i) = eigvec(i,i) + v(i)
  END DO
  CALL dsyev('v','l',n,eigvec,n,eig,work,5*n,info)
  DEALLOCATE(work) 
  srf(1:n) = eigvec(n,1:n) * srf_prm
  title='eigenvalues'
  CALL prntrm(title,eig,n,1,n,1,iout)
END SUBROUTINE Make_Hamiltonian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
