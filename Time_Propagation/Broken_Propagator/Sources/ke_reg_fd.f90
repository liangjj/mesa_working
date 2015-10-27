!deck ke_reg_fd.f
!***begin prologue     ke_reg_fd
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            fill regional finite differenc matrices from global ones
!***references
!***routines called
!***end prologue       ke__reg_fd
  SUBROUTINE ke_reg_fd(k_global,k_region,n_region,n_global,reg)
  USE prop_global
  USE Iterative_Global
  USE dvr_global,        ONLY: nreg
  USE dvrprop_global,    ONLY: title
  USE fd_global
  IMPLICIT NONE
  INTEGER                                :: n_region, n_global, reg
  REAL*8, DIMENSION(n_global,*)          :: k_global
  REAL*8, DIMENSION(n_region,n_region)   :: k_region
  INTEGER                                :: i
  CHARACTER (LEN=3)                      :: itoc
!
! Fill in the Kinetic Energy for all the special cases.
!
  DO i=1,n_region
     k_region(i,i) = k_global(1,i)
  END DO
  IF(n_region==2) THEN
     k_region(1,2) = dscale*d(2)
     k_region(2,1) = k_region(1,2)
  ELSE IF(n_region==3) THEN
     k_region(1,2) = dscale*d(2)
     k_region(2,1) = k_region(1,2)
     k_region(1,3) = dscale*d(3)
     k_region(3,1) = k_region(1,3)
     k_region(2,3) = k_region(1,2)
     k_region(3,2) = k_region(2,3)
  ELSE IF(n_region==4) THEN
     k_region(1,2) = dscale*d(2)
     k_region(2,1) = k_region(1,2)
     k_region(1,3) = dscale*d(3)
     k_region(3,1) = k_region(1,3)
     k_region(1,4) = dscale*d(4)
     k_region(4,1) = k_region(1,4)
     k_region(2,3) = k_region(1,2)
     k_region(3,2) = k_region(2,3)
     k_region(2,4) = k_region(1,3)
     k_region(4,2) = k_region(2,4)
     k_region(3,4) = k_region(1,2)
     k_region(4,3) = k_region(3,4)
  END IF
!
  k_region(1,1) = .5d0 * k_region(1,1)
  k_region(n_region,n_region) = .5d0 * k_region(n_region,n_region)
  IF(reg == 1) THEN
     k_region(1,1) = k_region(1,1)                                &
                          +                                       &
                     k_region(1,1)
  END If
  IF(reg == nreg) THEN
     k_region(n_region,n_region) = k_region(n_region,n_region)    &
                                              +                   &
                                   k_region(n_region,n_region)
  END IF
  IF (log_main(1)) THEN
      title= 'Modified Kinetic Energy Matrix Elements for '//  &
             ' Region = '//itoc(reg)
      call prntfmn(title,k_region,n_region,n_region,n_region,n_region,iout,'e')
  END IF
END SUBROUTINE ke_reg_fd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
