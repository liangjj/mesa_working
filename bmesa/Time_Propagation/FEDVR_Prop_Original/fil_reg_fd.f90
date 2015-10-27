!deck fil_reg_fd.f
!***begin prologue     fil_reg_fd
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            fill regional finite differenc matrices from global ones
!***references
!***routines called
!***end prologue       fil__reg_fd
  SUBROUTINE fil_reg_fd(kdiag,tr,nr,reg)
  USE arnoldi_global
  USE dvr_global,        ONLY: nreg
  USE dvrprop_global,    ONLY: title, keep_diag
  USE fd_global
  IMPLICIT NONE
  INTEGER                                :: nr, reg
  REAL*8, DIMENSION(nr)                  :: kdiag
  REAL*8, DIMENSION(nr,nr)               :: tr
  INTEGER                                :: i
  CHARACTER (LEN=3)                      :: itoc
!
! Fill in the Kinetic Energy for all the special cases.
!
  DO i=1,nr
     tr(i,i) = kdiag(i)
  END DO
  IF(nr==2) THEN
     tr(1,2) = dscale*d(2)
     tr(2,1) = tr(1,2)
  ELSE IF(nr==3) THEN
     tr(1,2) = dscale*d(2)
     tr(2,1) = tr(1,2)
     tr(1,3) = dscale*d(3)
     tr(3,1) = tr(1,3)
     tr(2,3) = tr(1,2)
     tr(3,2) = tr(2,3)
  ELSE IF(nr==4) THEN
     tr(1,2) = dscale*d(2)
     tr(2,1) = tr(1,2)
     tr(1,3) = dscale*d(3)
     tr(3,1) = tr(1,3)
     tr(1,4) = dscale*d(4)
     tr(4,1) = tr(1,4)
     tr(2,3) = tr(1,2)
     tr(3,2) = tr(2,3)
     tr(2,4) = tr(1,3)
     tr(4,2) = tr(2,4)
     tr(3,4) = tr(1,2)
     tr(4,3) = tr(3,4)
  END IF
!
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
END SUBROUTINE fil_reg_fd
