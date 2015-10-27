!
MODULE potential
!***begin prologue     potential
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           potential parameters
!***author             schneider, b. i.(nsf)
!***source             Modules
!***purpose            global shared variables for potential
!***description        data needed to describe time independent
!***                   and time dependent potentials
!
!***references

!***routines called    
!***end prologue       potential
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
  REAL*8                                  :: omega, width, shift, scale
  REAL*8                                  :: omega_t, width_t, shift_t, &
                                             scale_t
  REAL*8                                  :: gamma, delt, nl_coef
  REAL*8, DIMENSION(3)                    :: alpha, sigma, beta, x_0  
  REAL*8, DIMENSION(3)                    :: omega_q 
  REAL*8                                  :: d_0, e_0, t_0 
  REAL*8                                  :: natoms, scattering_length 
  REAL*8                                  :: nu_dw_plus, nu_dw_minus
  REAL*8                                  :: x_initial, x_final
  REAL*8                                  :: t_initial, t_final
  REAL*8                                  :: gamma_dw
  INTEGER                                 :: n_scale
  REAL*8                                  :: e_c, tau_c, k_0
  DATA alpha / 3*1.d0 /
  DATA sigma / 3*1.d0 /
  DATA x_0 / 3*0.d0 /
  DATA beta / 3*1.d0 /   
!
!
END MODULE potential
