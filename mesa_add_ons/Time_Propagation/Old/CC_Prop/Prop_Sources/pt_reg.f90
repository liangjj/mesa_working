!deck pt_reg.f
!***begin prologue     pt_reg
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           
!***author             schneider, barry (nsf)
!***source
!***purpose            fill regional points from global ones
!***references
!***routines called
!***end prologue       ke_reg_dvr
  SUBROUTINE pt_reg(pt_global,pt_region,n_region,n_global)
  USE Iterative_Global
  IMPLICIT NONE
  INTEGER                                :: n_region, n_global
  REAL*8, DIMENSION(n_global)            :: pt_global
  REAL*8, DIMENSION(n_region)            :: pt_region
!
  pt_region(1:n_region) =pt_global(1:n_region)
END SUBROUTINE pt_reg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
