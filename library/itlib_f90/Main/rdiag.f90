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
  SUBROUTINE rdiag
  USE io
  USE dvd_global
  USE dvd_prnt
  IMPLICIT NONE
  INTEGER                                :: info
  CHARACTER (LEN=3)                      :: itoc
  CALL dsyev('v','l',size,bwrk,maxvec,eigwrk,work,5*size,info)
  IF(info /= 0) THEN
     CALL lnkerr('error from direct diagonalization routine')
  END IF
  svec(1:size,1:size) = bwrk(1:size,1:size)
  IF(log_dvd(5)) THEN
     work(1:size,1) = eigwrk(1:size) + rep
     title='eigenvalues of small matrix iteration = ' //itoc(iter)
     CALL prntfm(title,work,size,1,size,1,iout)
  END IF
END SUBROUTINE rdiag

