!************************************************************************************
!deck read_hermite
!***begin prologue     read_hermite
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       read_hermite
  SUBROUTINE Read_Hermite
  USE dvr_global
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                               :: intkey
  LOGICAL                               :: logkey
  nreg=1
  nfix=0
  fix(1)=.false.
  fix(2)=.false.
  drop(1)=.false.
  drop(2)=.false.
  CALL intarr(card,'number-of-points',n,nreg,' ')
  npt(1)=n(1)
!*************************************************************************************
  END SUBROUTINE Read_Hermite
