!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            MODULE so_exponentiation
                            USE dvrprop_global
                            USE dvr_shared
                            USE so_exponential_diagonal_propagator
                            USE so_exponential_diagonal_multiplication
                            USE so_exponential_off_diagonal_multiplication
                            USE plot_wavefunction
                            USE plot_propagator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_exponentiation
!***begin prologue     so_exponentiation
!***date written       040707   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            solve time-dependent schroedinger equation using a
!***                   second order accurate split-operator method with
!***                   a DVR or FD representation.
!***
!***description        driver for the calculation of the time propagation.
!                      the diagonal exponential scaling is performed in this
!                      subroutine but the real work occurs in exp_off_diag.
!***references
!***routines called    
!***end prologue       so_exponentiation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        INTERFACE so_exp
                MODULE PROCEDURE so_exp_d,                               &
                                 so_exp_z
                    END INTERFACE so_exp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_exp_d.f
!***begin prologue     so_exp_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            time-dependent wave function calculation
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       so_exp_d
  SUBROUTINE so_exp_d(wave_function,scratch_vector)
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)       :: wave_function
  REAL*8, DIMENSION(n3d)       :: scratch_vector
!
! Since the potential is, in general, time-dependent, diagonal
! propagators will be reconstructed at each time-step.
!
!
  scratch_vector(:) = 0.d0
  IF(log_main(9)) then
     title='Initial vector'
     CALL print_psi(wave_function)
  END IF
  IF (prop_order  == 2 ) THEN
!
!            First Diagonal Scaling by Potential
!                t = p_fac * deltat *.5 / hbar
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL so_exp_diagonal_propagator(exp_diag_d)
      IF(log_main(9)) then
          title='Exponential Diagonal Propagator'
         CALL plot_prop(exp_diag_d)
      END IF
!
      CALL so_exp_diagonal_multiply(wave_function)
!
     IF(log_main(9)) then
         title='First Diagonally Scaled Vector'
         CALL print_psi(wave_function)
      END IF
!
!            Off diagonal scaling
!                t = p_fac * deltat / hbar
      prop_point = 1
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
      IF(log_main(9)) then
         title='Off Diagonally Propagated Vector'
         CALL print_psi(wave_function)
      END IF
!
!            Second Diagonal Scaling by Potential
!                t = p_fac * deltat *.5 / hbar
!
      CALL so_exp_diagonal_multiply(wave_function)
!
      IF(log_main(9)) then
         title='Second Diagonally Scaled Vector'
         CALL print_psi(wave_function)
      END IF
  ELSE IF (prop_order == 4) THEN
!
!            First Diagonal Scaling by Potential
!               t = p_fac * deltat * .5 / h_bar  
!
      tau_loc = p_fac * deltat * .5d0 / hbar
!
      CALL so_exp_diagonal_propagator(exp_diag_d)
      CALL so_exp_diagonal_multiply(wave_function)
!
!            First Off-Diagonal Scaling 
!            The second order split operator is used
!               t = p_fac * deltat / hbar
!
      prop_point = 1
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Second Diagonal Scaling by potential
!               t = p_fac * deltat / hbar 
!
      tau_loc = p_fac * deltat / hbar
      CALL so_exp_diagonal_propagator(exp_diag_d)
      CALL so_exp_diagonal_multiply(wave_function)
!
!            Second Off-Diagonal Scaling
!            Same as first
!
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Third Diagonal Scaling
!               t = ( 1. - 3. * p_fac ) * deltat *.5 / hbar
!
      tau_loc = ( 1.d0 - 3.d0 * p_fac ) * deltat * .5d0 / hbar
      CALL so_exp_diagonal_propagator(exp_diag_d)
      CALL so_exp_diagonal_multiply(wave_function)
!
!            Third Off_diagonal Scaling
!               t = ( 1. - 4. * p_fac ) * deltat / hbar    
!
      prop_point = 2
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Fourth Diagonal Scaling
!               Same as the third.
!
      CALL so_exp_diagonal_multiply(wave_function)
!
!            Fourth Off-Diagonal Scaling
!              Same  as first.
!
      prop_point = 1
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Fifth Diagonal Scaling
!            Same as second
!
      tau_loc = p_fac * deltat / hbar
      CALL so_exp_diagonal_propagator(exp_diag_d)
      CALL so_exp_diagonal_multiply(wave_function)
!
!            Fifth Off-Diagonal Scaling
!            Same as first
!
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Sixth Diagonal Scaling
!            Same as first
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL so_exp_diagonal_propagator(exp_diag_d)
      CALL so_exp_diagonal_multiply(wave_function)
!
!     OK.  We are done.
! 
  END IF
END SUBROUTINE so_exp_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck so_exp_z.f
!***begin prologue     so_exp_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            time-dependent wave function calculation
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       so_exp_3d_z
  SUBROUTINE so_exp_z(wave_function,scratch_vector)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d)       :: wave_function
  COMPLEX*16, DIMENSION(n3d)       :: scratch_vector
!
! Since the potential is, in general, time-dependent, diagonal
! propagators will be reconstructed at each time-step.
!
!
  scratch_vector(:) = (0.d0,0.d0)
  IF(log_main(9)) then
     title='Initial vector'
     CALL print_psi(wave_function)
  END IF
  IF (prop_order  == 2 ) THEN
!
!            First Diagonal Scaling by Potential
!                t = p_fac * deltat *.5 / hbar
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL so_exp_diagonal_propagator(exp_diag_z)
      IF(log_main(9)) then
          title='Exponential Diagonal Propagator'
         CALL plot_prop(exp_diag_z)
      END IF
!
      CALL so_exp_diagonal_multiply(wave_function)
!
     IF(log_main(9)) then
        title='First Diagonally Scaled Vector'
        CALL print_psi(wave_function)
     END IF
!
!            Off diagonal scaling
!                t = p_fac * deltat / hbar
     prop_point = 1
     CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
     IF(log_main(9)) then
        title='Off Diagonally Propagated Vector'
        CALL print_psi(wave_function)
     END IF
!
!            Second Diagonal Scaling by Potential
!                t = p_fac * deltat *.5 / hbar
!
     CALL so_exp_diagonal_multiply(wave_function)
!
     IF(log_main(9)) then
        title='Second Diagonally Scaled Vector'
        CALL print_psi(wave_function)
     END IF
  ELSE IF (prop_order == 4) THEN
!
!            First Diagonal Scaling by Potential
!               t = p_fac * deltat * .5 / h_bar  
!
      tau_loc = p_fac * deltat * .5d0 / hbar
!
      CALL so_exp_diagonal_propagator(exp_diag_z)
      CALL so_exp_diagonal_multiply(wave_function)
!
!            First Off-Diagonal Scaling 
!            The second order split operator is used
!               t = p_fac * deltat / hbar
!
      prop_point = 1
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Second Diagonal Scaling by potential
!               t = p_fac * deltat / hbar 
!
      tau_loc = p_fac * deltat / hbar
      CALL so_exp_diagonal_propagator(exp_diag_z)
      CALL so_exp_diagonal_multiply(wave_function)
!
!            Second Off-Diagonal Scaling
!            Same as first
!
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Third Diagonal Scaling
!               t = ( 1. - 3. * p_fac ) * deltat *.5 / hbar
!
      tau_loc = ( 1.d0 - 3.d0 * p_fac ) * deltat * .5d0 / hbar
      CALL so_exp_diagonal_propagator(exp_diag_z)
      CALL so_exp_diagonal_multiply(wave_function)
!
!            Third Off_diagonal Scaling
!               t = ( 1. - 4. * p_fac ) * deltat / hbar    
!
      prop_point = 2
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Fourth Diagonal Scaling
!               Same as the third.
!
      CALL so_exp_diagonal_multiply(wave_function)
!
!            Fourth Off-Diagonal Scaling
!              Same  as first.
!
      prop_point = 1
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Fifth Diagonal Scaling
!            Same as second
!
      tau_loc = p_fac * deltat / hbar
      CALL so_exp_diagonal_propagator(exp_diag_z)
      CALL so_exp_diagonal_multiply(wave_function)
!
!            Fifth Off-Diagonal Scaling
!            Same as first
!
      CALL so_exp_off_diagonal_multiply(wave_function,scratch_vector)
!
!            Sixth Diagonal Scaling
!            Same as first
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL so_exp_diagonal_propagator(exp_diag_z)
      CALL so_exp_diagonal_multiply(wave_function)
!
!     OK.  We are done.
! 
  END IF
END SUBROUTINE so_exp_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE so_exponentiation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
