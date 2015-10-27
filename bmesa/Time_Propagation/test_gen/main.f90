  program main
  USE prop_umat_2
  USE prop_diag_2
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  INTEGER                                :: m, reg
  REAL*8, DIMENSION(10)                  :: eval
  REAL*8, DIMENSION(10,10)               :: k_mat
  COMPLEX*16, DIMENSION(10,10)           :: k_mat_c
  REAL*8, DIMENSION(10,10)               :: evec_r
  REAL*8, DIMENSION(10,10)               :: evec_l
  REAL*8, DIMENSION(10,10)               :: cos_t, sin_t
  COMPLEX*16, DIMENSION(10)              :: eval_c
  COMPLEX*16, DIMENSION(10,10)           :: evec_c_r
  COMPLEX*16, DIMENSION(10,10)           :: evec_c_l
  REAL*8, DIMENSION(10,10)               :: si, ci
  COMPLEX*16, DIMENSION(10,10)           :: si_c, ci_c
  REAL*8                                 :: tau
  m=10
  reg=1
  call diag_reg(k_mat,eval,evec_r,evec_l,m,reg)
  call diag_reg(k_mat_c,eval_c,evec_c_r,evec_c_l,m,reg)
  call umat_reg(eval,evec_r,evec_l,cos_t,sin_t,si,ci,tau,m)
  call umat_reg(eval_c,evec_c_r,evec_c_l,cos_t,sin_t,  &
                si_c,ci_c,tau,m)
  stop
END PROGRAM main
