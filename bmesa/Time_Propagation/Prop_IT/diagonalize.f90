!***********************************************************************
                           MODULE diagonalize
                           INTERFACE h_diag
                    MODULE PROCEDURE h_diag_d, h_diag_z
                       END INTERFACE h_diag
!***********************************************************************
                           CONTAINS
!***********************************************************************
!***********************************************************************
!***********************************************************************
!deck h_diag_z
!**begin prologue     h_diag_z
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           diagonalization
!**author             schneider, barry (nsf)
!**source
!**purpose            Driver for hermitian diagonalization.
!**references
!**routines called
!**end prologue       h_diag_z
  SUBROUTINE h_diag_z(small_matrix_s)
  USE arnoldi_global_rt
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(maxvec,maxvec)   :: small_matrix_s
  INTEGER                                :: info
  CHARACTER (LEN=3)                      :: itoc

! Diagonalize the Hermitian matrix for the real eigenvalues.

  CALL zheev('v','l',size,small_matrix_s,maxvec,eig,work,     &
                   5*size,rwork,info)
  IF(info /= 0) THEN
     CALL lnkerr('error from direct diagonalization routine')
  END IF
  title='eigenvalues of small matrix iteration = ' //itoc(iter)
  CALL prntfm(title,eig,size,1,size,1,iout)
  IF(log_prp(6)) THEN
       title='eigenvectors of small matrix iteration = ' //itoc(iter)
       CALL prntcm(title,small_matrix_s,size,size,maxvec,maxvec,iout)
  END IF
END SUBROUTINE h_diag_z
!deck h_diag_d
!**begin prologue     h_diag_d
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           diagonalization
!**author             schneider, barry (nsf)
!**source
!**purpose            driver for real, symmetric diagonalization.
!**references
!**routines called
!**end prologue       h_diag_d
  SUBROUTINE h_diag_d(small_matrix_s)
  USE arnoldi_global_it
  IMPLICIT NONE
  REAL*8                                 :: eig_tst
  REAL*8, DIMENSION(maxvec,maxvec)       :: small_matrix_s
  INTEGER                                :: info
  CHARACTER (LEN=3)                      :: itoc

! Diagonalize the Hermitian matrix for the real eigenvalues.

  CALL dsyev('v','l',size,small_matrix_s,maxvec,eig,work,     &
                   5*size,rwork,info)
  IF(info /= 0) THEN
     CALL lnkerr('error from direct diagonalization routine')
  END IF
  title='eigenvalues of small matrix iteration = ' //itoc(iter)
  CALL prntfm(title,eig,size,1,size,1,iout)
  IF(log_prp(6)) THEN
       title='eigenvectors of small matrix iteration = ' //itoc(iter)
       CALL prntrm(title,small_matrix_s,size,size,maxvec,maxvec,iout)
  END IF
END SUBROUTINE h_diag_d
END MODULE diagonalize

