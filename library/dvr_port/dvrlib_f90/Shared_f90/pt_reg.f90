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
!***end prologue       pt_reg_dvr
  SUBROUTINE pt_reg(pt_global,wt_global,pt_region,nrm_region,n_region,n_global)
  IMPLICIT NONE
  INTEGER                                :: n_region, n_global
  REAL*8, DIMENSION(n_global)            :: pt_global, wt_global
  REAL*8, DIMENSION(n_region)            :: pt_region, nrm_region
  INTEGER                                :: i
  CHARACTER (LEN=3)                      :: itoc
!
  pt_region(1:n_region) = pt_global(1:n_region)
  nrm_region(1:n_region) = 1.d0/sqrt(wt_global(1:n_region))
END SUBROUTINE pt_reg
