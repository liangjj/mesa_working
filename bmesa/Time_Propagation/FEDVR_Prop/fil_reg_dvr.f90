!deck fil_reg_dvr.f
!***begin prologue     fil_reg_dvr
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            fill regional DVR matrices from global ones
!***references
!***routines called
!***end prologue       fil_reg_dvr
  SUBROUTINE fil_reg_dvr(kmat,tr,nr,n,reg)
  USE arnoldi_global
  USE dvr_global,        ONLY: nreg
  USE dvrprop_global,    ONLY: keep_diag
  IMPLICIT NONE
  INTEGER                                :: nr, n, reg
  REAL*8, DIMENSION(n,*)                 :: kmat
  REAL*8, DIMENSION(nr,nr)               :: tr
  INTEGER                                :: i
  CHARACTER (LEN=3)                      :: itoc
!
  tr =kmat(1:nr,1:nr)
  tr(1,1) = .5d0 * tr(1,1)
  tr(nr,nr) = .5d0 * tr(nr,nr)
  IF(reg == 1) THEN
     tr(1,1) = tr(1,1) + tr(1,1)
  END If
  IF(reg == nreg) THEN
     tr(nr,nr) = tr(nr,nr) + tr(nr,nr)
  END IF
  IF (log_main(1)) THEN
      title= 'Modified Kinetic Energy Matrix Elements for '//  &
             ' Region = '//itoc(reg)
      call prntfmn(title,tr,nr,nr,nr,nr,iout,'e')
  END IF
END SUBROUTINE fil_reg_dvr
