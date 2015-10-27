!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE regional_diag
!**begin prologue     regional_diag_2
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       regional_diag_2
!
  INTERFACE diag_reg
  SUBROUTINE diag_reg_d(k_region,eval_region,evec_region_right,  &
                        evec_region_left,n_region,reg)
  USE dvr_global
  USE dvrprop_global
  IMPLICIT NONE
  INTEGER                                :: n_region, reg
  INTEGER                                :: info
  REAL*8, DIMENSION(n_region,n_region)   :: k_region
  REAL*8, DIMENSION(n_region)            :: eval_region
  REAL*8, DIMENSION(n_region,n_region)   :: evec_region_right
  REAL*8, DIMENSION(n_region,n_region)   :: evec_region_left
  CHARACTER(LEN=4)                       :: itoc
  END SUBROUTINE diag_reg_d  
  SUBROUTINE diag_reg_z(k_region,eval_region,evec_region_right, &
                        evec_region_left,n_region,reg)
  USE dvr_global
  USE dvrprop_global
  IMPLICIT NONE
  INTEGER                                  :: n_region, reg
  INTEGER                                  :: info
  INTEGER                                  :: i, j
  COMPLEX*16, DIMENSION(n_region,n_region) :: k_region
  COMPLEX*16, DIMENSION(n_region)          :: eval_region
  COMPLEX*16, DIMENSION(n_region,n_region) :: evec_region_right
  COMPLEX*16, DIMENSION(n_region,n_region) :: evec_region_left
  COMPLEX*16                               :: val_z, cdotc
END SUBROUTINE diag_reg_z
END INTERFACE diag_reg
END MODULE regional_diag

