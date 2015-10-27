!deck rdiag.f
!***begin prologue     rdiag
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           diagonalization
!***author             schneider, barry (nsf)
!***source
!***purpose            driver for real diagonalization.
!***
!***references

!***routines called
!**end prologue       rdiag
  SUBROUTINE rdiag(rep,iter,size)
  USE io
  USE dvd_global
  USE dvd_prnt
  IMPLICIT NONE
  REAL*8                                 :: rep
  INTEGER                                :: iter
  INTEGER                                :: size
  INTEGER                                :: info
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  CALL dsyev('v','l',size,mat,maxvec,eig,work,5*size,info)
  IF(info /= 0) THEN
     CALL lnkerr('error from direct diagonalization routine')
  END IF
  vec(1:size,1:size) = mat(1:size,1:size)
  IF(log_dvd(5)) THEN
     work = eig + rep
     title='eigenvalues of small matrix iteration = ' //itoc(iter)
     CALL prntfm(title,work,size,1,size,1,iout)
  END IF
END SUBROUTINE rdiag

