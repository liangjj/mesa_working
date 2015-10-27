!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE regional_umat
!**begin prologue     regional_umat_2
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       regional_umat_2
!
  INTERFACE umat_reg
  SUBROUTINE umat_reg_d(eval_region,evec_region_right,    &
                        evec_region_left,cos_t,sin_t,     &
                        si,ci,tau,n_region)
  USE dvrprop_global_rt
  IMPLICIT NONE
  INTEGER                                     :: n_region
  INTEGER                                     :: i
  REAL*8, DIMENSION(n_region)                 :: eval_region
  REAL*8, DIMENSION(n_region,n_region)        :: evec_region_right
  REAL*8, DIMENSION(n_region,n_region)        :: evec_region_left
  REAL*8, DIMENSION(n_region,n_region)        :: cos_t, sin_t
  REAL*8, DIMENSION(n_region,n_region)        :: si, ci
  REAL*8                                      :: ex_fac, tau
  CHARACTER*16                                :: fptoc
  END SUBROUTINE umat_reg_d
  SUBROUTINE umat_reg_z(eval_region,evec_region_right,     &
                        evec_region_left,cos_t,sin_t,      &
                        si,ci,tau,n_region)
  USE dvrprop_global_rt
  IMPLICIT NONE
  INTEGER                                     :: n_region
  INTEGER                                     :: i
  COMPLEX*16, DIMENSION(n_region)             :: eval_region
  COMPLEX*16, DIMENSION(n_region,n_region)    :: evec_region_left
  COMPLEX*16, DIMENSION(n_region,n_region)    :: evec_region_right
  REAL*8,     DIMENSION(n_region,n_region)    :: cos_t, sin_t
  COMPLEX*16, DIMENSION(n_region,n_region)    :: si, ci
  REAL*8                                      :: tau, ex_fac
  REAL*8                                      :: re_eig, im_eig
  CHARACTER*16                                :: fptoc
  END SUBROUTINE umat_reg_z
  END INTERFACE umat_reg
END MODULE regional_umat
