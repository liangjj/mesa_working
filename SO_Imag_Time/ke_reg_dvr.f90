!deck ke_reg_dvr.f
!***begin prologue     ke_reg_dvr
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            fill regional DVR matrices from global ones
!***references
!***routines called
!***end prologue       ke_reg_dvr
  SUBROUTINE ke_reg_dvr(k_global,k_region,n_region,n_global,reg)
  USE arnoldi_global
  USE dvr_global,           ONLY: nreg
  IMPLICIT NONE
  INTEGER                                :: n_region, n_global, reg
  REAL*8, DIMENSION(n_global,*)          :: k_global
  REAL*8, DIMENSION(n_region,n_region)   :: k_region
  INTEGER                                :: i
  CHARACTER (LEN=3)                      :: itoc
!
  k_region =k_global(1:n_region,1:n_region)
  k_region(1,1) = .5d0 * k_region(1,1)
  k_region(n_region,n_region) = .5d0 * k_region(n_region,n_region)
  IF(reg == 1) THEN
     k_region(1,1) = k_region(1,1)                               &
                                   +                             &
                     k_region(1,1)
  END If
  IF(reg == nreg) THEN
     k_region(n_region,n_region) = k_region(n_region,n_region)   &
                                              +                  &
                                   k_region(n_region,n_region)
  END IF
  IF (log_main(1)) THEN
      title= 'Modified Kinetic Energy Matrix Elements for '//  &
             ' Region = '//itoc(reg)
      call prntfmn(title,k_region,n_region,n_region,n_region,n_region,iout,'e')
  END IF
END SUBROUTINE ke_reg_dvr
