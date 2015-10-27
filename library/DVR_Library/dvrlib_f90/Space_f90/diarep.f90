!deck diarep.f
!***begin prologue     diarep
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            diagonalization of DVR/FEM hamiltonian.
!***
!***references
!***routines called
!***end prologue       diarep
  SUBROUTINE diarep(h,eigv,eigvec,nphy)
  USE input_output
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nphy, info
  INTEGER                                :: i
  REAL*8, DIMENSION(nphy,nphy)           :: h, eigvec
  REAL*8, DIMENSION(nphy)                :: eigv
  REAL*8, DIMENSION(:), ALLOCATABLE      :: work
  CHARACTER (LEN=80)                     :: title
  eigvec=0.d0
  CALL copy(h,eigvec,nphy*nphy)
  IF(prn(9)) THEN
     title='physical hamiltonian'
     call prntrm(title,h,nphy,nphy,nphy,nphy,iout)
  END IF
  ALLOCATE(work(5*nphy))
  CALL dsyev('v','l',nphy,h,nphy,eigv,work,5*nphy,info)
  DEALLOCATE(work)
  IF(prn(9)) THEN
     title='eigenvalues'
     CALL prntrm(title,eigv,nphy,1,nphy,1,iout)
  END IF
  IF(prn(10)) THEN
     title='eigenvectors'
     CALL prntrm(title,eigvec,nphy,nphy,nphy,nphy,iout)
  END IF
END SUBROUTINE diarep



