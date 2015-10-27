!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        MODULE real_space_propagator
!deck real_space_propagator
!***begin prologue     real_space_propagator
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
!***routines called    exp_off_diag, prntrm
!***end prologue       real_space_propagator
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        INTERFACE rsp
                MODULE PROCEDURE rsp_1d_z, rsp_2d_z, rsp_3d_z
                    END INTERFACE rsp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck rsp_1d_z.f
!***begin prologue     rsp_1d_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           propagation
!***author             schneider, b. i.(nsf)
!***source
!***purpose            time-dependent wave function calculation
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       rsp_1d_z
  SUBROUTINE rsp_1d_z(wave_function,scratch_vector)
     USE dvrprop_global_rt
     USE dvr_shared
     USE exp_off_diagonal
     USE plot_wavefunction
     USE exp_off_diagonal
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)                       :: wave_function
  REAL*8, DIMENSION(nphy(1),2)                       :: scratch_vector

!
! Since the potential is, in general, time-dependent, diagonal
! propagators will be reconstructed at each time-step.
!
!
  IF(log_main(9)) then
     title='Initial vector'
     call print_psi(wave_function)
  END IF
  IF (prop_order  == 2 ) THEN
!----------------------------------------------------------------------
!
!            First Diagonal Scaling by Potential
!            t = p_fac * deltat *.5 / hbar
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL cs_diag
      IF(log_main(9)) then
         title='Diagonal Cosine Propagator'
         CALL plot_propagator(cos_diag)
         title='Diagonal Sine Propagator'
         CALL plot_propagator(sin_diag)
      END IF
!
      call diagonal_mul(wave_function,scratch_vector)
!
     IF(log_main(9)) then
        title='First Diagonally Scaled Vector'
        CALL print_psi(wave_function)
     END IF
!----------------------------------------------------------------------
!
!            Off diagonal scaling
!            t = p_fac * deltat / hbar
     prop_point = 1
     CALL exp_off_diag(wave_function,scratch_vector)
!
     IF(log_main(9)) then
        title='Off Diagonally Propagated Vector'
        CALL print_psi(wave_function)
     END IF
!
!----------------------------------------------------------------------
!            Second Diagonal Scaling by Potential
!            t = p_fac * deltat *.5 / hbar
!
     call diagonal_mul(wave_function,scratch_vector)
!
     IF(log_main(9)) then
        title='Second Diagonally Scaled Vector'
        CALL print_psi(wave_function)
     END IF
!
!----------------------------------------------------------------------
  ELSE IF (prop_order == 4) THEN
!
!----------------------------------------------------------------------
!            First Diagonal Scaling by Potential
!            t = p_fac * deltat * .5 / h_bar  
!
      tau_loc = p_fac * deltat * .5d0 / hbar
!
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!----------------------------------------------------------------------
!            First Off-Diagonal Scaling 
!            The second order split operator is used
!            t = p_fac * deltat / hbar
!
      prop_point = 1
      CALL exp_off_diag(wave_function,scratch_vector)
!
!----------------------------------------------------------------------
!            Second Diagonal Scaling by potential
!            t = p_fac * deltat / hbar 
!
      tau_loc = p_fac * deltat / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!----------------------------------------------------------------------
!            Second Off-Diagonal Scaling
!            Same as first
!
      CALL exp_off_diag(wave_function,scratch_vector)
!
!----------------------------------------------------------------------
!            Third Diagonal Scaling
!            t = ( 1. - 3. * p_fac ) * deltat *.5 / hbar
!
      tau_loc = ( 1.d0 - 3.d0 * p_fac ) * deltat * .5d0 / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!----------------------------------------------------------------------
!            Third Off_diagonal Scaling
!               t = ( 1. - 4. * p_fac ) * deltat / hbar    
!
      prop_point = 2
      CALL exp_off_diag(wave_function,scratch_vector)
!
!----------------------------------------------------------------------
!            Fourth Diagonal Scaling
!               Same as the third.
!
      call diagonal_mul(wave_function,scratch_vector)
!
!----------------------------------------------------------------------
!            Fourth Off-Diagonal Scaling
!              Same  as first.
!
      prop_point = 1
      CALL exp_off_diag(wave_function,scratch_vector)
!
!----------------------------------------------------------------------
!            Fifth Diagonal Scaling
!            Same as second
!
      tau_loc = p_fac * deltat / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!----------------------------------------------------------------------
!            Fifth Off-Diagonal Scaling
!            Same as first
!
      CALL exp_off_diag(wave_function,scratch_vector)
!
!----------------------------------------------------------------------
!            Sixth Diagonal Scaling
!            Same as first
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!----------------------------------------------------------------------
!                  OK.  We are done.
! 
  END IF
END SUBROUTINE rsp_1d_z
!deck rsp_2d_z.f
!***begin prologue     rsp_2d_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            time-dependent wave function calculation
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       rsp_2d_z
  SUBROUTINE rsp_2d_z(wave_function,scratch_vector)
     USE dvrprop_global_rt
     USE dvr_shared
     USE exp_off_diagonal
     USE plot_wavefunction
     USE exp_off_diagonal
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),2)                  :: wave_function
  REAL*8, DIMENSION(nphy(2),nphy(1),2)                  :: scratch_vector
!
! Since the potential is, in general, time-dependent, diagonal
! propagators will be reconstructed at each time-step.
!
!
  IF(log_main(9)) then
     title='Initial vector'
     call print_psi(wave_function)
  END IF
  IF (prop_order  == 2 ) THEN
!
!            First Diagonal Scaling by Potential
!                t = p_fac * deltat *.5 / hbar
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL cs_diag
      IF(log_main(9)) then
         title='Diagonal Cosine Propagator'
         CALL plot_propagator(cos_diag)
         title='Diagonal Sine Propagator'
         CALL plot_propagator(sin_diag)
      END IF
!
      call diagonal_mul(wave_function,scratch_vector)
!
     IF(log_main(9)) then
        title='First Diagonally Scaled Vector'
        CALL print_psi(wave_function)
     END IF
!
!            Off diagonal scaling
!                t = p_fac * deltat / hbar
     prop_point = 1
     CALL exp_off_diag(wave_function,scratch_vector)
!
     IF(log_main(9)) then
        title='Off Diagonally Propagated Vector'
        CALL print_psi(wave_function)
     END IF
!
!            Second Diagonal Scaling by Potential
!                t = p_fac * deltat *.5 / hbar
!
     call diagonal_mul(wave_function,scratch_vector)
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
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!            First Off-Diagonal Scaling 
!            The second order split operator is used
!               t = p_fac * deltat / hbar
!
      prop_point = 1
      CALL exp_off_diag(wave_function,scratch_vector)
!
!            Second Diagonal Scaling by potential
!               t = p_fac * deltat / hbar 
!
      tau_loc = p_fac * deltat / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!            Second Off-Diagonal Scaling
!            Same as first
!
      CALL exp_off_diag(wave_function,scratch_vector)
!
!            Third Diagonal Scaling
!               t = ( 1. - 3. * p_fac ) * deltat *.5 / hbar
!
      tau_loc = ( 1.d0 - 3.d0 * p_fac ) * deltat * .5d0 / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!            Third Off_diagonal Scaling
!               t = ( 1. - 4. * p_fac ) * deltat / hbar    
!
      prop_point = 2
      CALL exp_off_diag(wave_function,scratch_vector)
!
!            Fourth Diagonal Scaling
!               Same as the third.
!
      call diagonal_mul(wave_function,scratch_vector)
!
!            Fourth Off-Diagonal Scaling
!              Same  as first.
!
      prop_point = 1
      CALL exp_off_diag(wave_function,scratch_vector)
!
!            Fifth Diagonal Scaling
!            Same as second
!
      tau_loc = p_fac * deltat / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!            Fifth Off-Diagonal Scaling
!            Same as first
!
      CALL exp_off_diag(wave_function,scratch_vector)
!
!            Sixth Diagonal Scaling
!            Same as first
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!     OK.  We are done.
! 
  END IF

END SUBROUTINE rsp_2d_z
!deck rsp_3d_z.f
!***begin prologue     rsp_3d_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            time-dependent wave function calculation
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       rsp_3d_z
  SUBROUTINE rsp_3d_z(wave_function,scratch_vector)
     USE dvrprop_global_rt
     USE dvr_shared
     USE exp_off_diagonal
     USE plot_wavefunction
     USE exp_off_diagonal
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)       :: wave_function
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)       :: scratch_vector
!
! Since the potential is, in general, time-dependent, diagonal
! propagators will be reconstructed at each time-step.
!
!
  IF(log_main(9)) then
     title='Initial vector'
     call print_psi(wave_function)
  END IF
  IF (prop_order  == 2 ) THEN
!
!            First Diagonal Scaling by Potential
!                t = p_fac * deltat *.5 / hbar
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL cs_diag
      IF(log_main(9)) then
         title='Diagonal Cosine Propagator'
         CALL plot_propagator(cos_diag)
         title='Diagonal Sine Propagator'
         CALL plot_propagator(sin_diag)
      END IF
!
      call diagonal_mul(wave_function,scratch_vector)
!
     IF(log_main(9)) then
        title='First Diagonally Scaled Vector'
        CALL print_psi(wave_function)
     END IF
!
!            Off diagonal scaling
!                t = p_fac * deltat / hbar
     prop_point = 1
     CALL exp_off_diag(wave_function,scratch_vector)
!
     IF(log_main(9)) then
        title='Off Diagonally Propagated Vector'
        CALL print_psi(wave_function)
     END IF
!
!            Second Diagonal Scaling by Potential
!                t = p_fac * deltat *.5 / hbar
!
     call diagonal_mul(wave_function,scratch_vector)
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
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!            First Off-Diagonal Scaling 
!            The second order split operator is used
!               t = p_fac * deltat / hbar
!
      prop_point = 1
      CALL exp_off_diag(wave_function,scratch_vector)
!
!            Second Diagonal Scaling by potential
!               t = p_fac * deltat / hbar 
!
      tau_loc = p_fac * deltat / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!            Second Off-Diagonal Scaling
!            Same as first
!
      CALL exp_off_diag(wave_function,scratch_vector)
!
!            Third Diagonal Scaling
!               t = ( 1. - 3. * p_fac ) * deltat *.5 / hbar
!
      tau_loc = ( 1.d0 - 3.d0 * p_fac ) * deltat * .5d0 / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!            Third Off_diagonal Scaling
!               t = ( 1. - 4. * p_fac ) * deltat / hbar    
!
      prop_point = 2
      CALL exp_off_diag(wave_function,scratch_vector)
!
!            Fourth Diagonal Scaling
!               Same as the third.
!
      call diagonal_mul(wave_function,scratch_vector)
!
!            Fourth Off-Diagonal Scaling
!              Same  as first.
!
      prop_point = 1
      CALL exp_off_diag(wave_function,scratch_vector)
!
!            Fifth Diagonal Scaling
!            Same as second
!
      tau_loc = p_fac * deltat / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!            Fifth Off-Diagonal Scaling
!            Same as first
!
      CALL exp_off_diag(wave_function,scratch_vector)
!
!            Sixth Diagonal Scaling
!            Same as first
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!     OK.  We are done.
! 
  END IF
END SUBROUTINE rsp_3d_z
END MODULE real_space_propagator
