!************************************************************************************
!deck read_laguerre
!***begin prologue     read_laguerre
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           dvr, basis, orthogonal polynomial, input
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            input subroutine for dvr basis sets
!***description        user interface for dvr library
!***references
!***routines called    
!***end prologue       read_laguerre
  SUBROUTINE Read_Laguerre
  USE dvr_global
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                               :: intkey
  LOGICAL                               :: logkey
  nreg=1
  nfix=intkey(card,'number_of_fixed_points',0,' ')
  CALL intarr(card,'number_of_points',n,nreg,' ')
  IF(nfix == 1) THEN
     fix(1)=logkey(card,'left_fixed_point',.true.,' ')
     drop(1)=logkey(card,'drop_left_function',.false.,' ')
  END IF
  bcl=1
  IF(drop(1)) THEN
     bcl=0
  END IF
  npt(1)=n(1)
!*************************************************************************************
  END SUBROUTINE Read_Laguerre
