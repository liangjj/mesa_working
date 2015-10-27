!deck tri_diag
!**begin prologue     tri_diag
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           diagonalization
!**author             schneider, barry (nsf)
!**source
!**purpose            eigenvalues and eigenvectors of tridiagonal matrix.
!**references
!**routines called
!**end prologue       tri_diag
  SUBROUTINE tri_diag(it)
  USE lanczos_global
  IMPLICIT NONE
  INTEGER                                :: it, jt, ierr
  CHARACTER (LEN=3)                      :: itoc

! Diagonalize the tridiagonal matrix for the real eigenvalues.
  b(1:it,1:it) = 0.d0
  DO jt=1,it
     b(jt,jt) = 1.d0
  END DO
  eig(1:it) = a_lanczos(1:it)
  bwrk(1:it) = b_lanczos(1:it)
!  IF(log_prp(5)) THEN
!     title='lanczos a coefficients iteration = ' //itoc(it)
!     CALL prntrm(title,a_lanczos,it,1,it,1,iout)
!     title='lanczos b coefficients iteration = ' //itoc(it)
!     CALL prntrm(title,b_lanczos,it,1,it,1,iout)
!  END IF
  call  imtql2(maxit,it,eig,bwrk,b,ierr)
  IF(ierr /= 0) THEN
     CALL lnkerr('error from direct diagonalization routine')
  END IF
  IF(log_prp(5)) THEN
     title='eigenvalues of small matrix iteration = ' //itoc(it)
     CALL prntfm(title,eig,it,1,it,1,iout)
     title='eigenvectors of small matrix iteration = ' //itoc(it)
     CALL prntrm(title,b,it,it,maxit,maxit,iout)
  END IF
END SUBROUTINE tri_diag

