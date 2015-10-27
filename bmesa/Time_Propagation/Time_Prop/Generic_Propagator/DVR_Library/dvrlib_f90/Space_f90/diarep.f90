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
  SUBROUTINE diarep(h,eigv,eigvec,work,nphy)
  USE input_output
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nphy, info
  REAL*8, DIMENSION(nphy,nphy)           :: h, eigvec
  REAL*8, DIMENSION(nphy)                :: eigv
  REAL*8, DIMENSION(5*nphy)              :: work
  CHARACTER (LEN=80)                     :: title
  CALL copy(h,eigvec,nphy*nphy)
  IF(prn(9)) THEN
     title='physical hamiltonian'
     call prntrm(title,h,nphy,nphy,nphy,nphy,iout)
  END IF
  CALL dsyev('v','l',nphy,eigvec,nphy,eigv,work,5*nphy,info)
  IF(prn(9)) THEN
     title='eigenvalues'
     CALL prntrm(title,eigv,nphy,1,nphy,1,iout)
  END IF
  IF(prn(10)) THEN
     title='eigenvectors'
     CALL prntrm(title,eigvec,nphy,nphy,nphy,nphy,iout)
  END IF
END SUBROUTINE diarep



