!deck rsp_1d.f
!***begin prologue     rsp_1d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            time-dependent wave function calculation
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       rsp_1d
  SUBROUTINE rsp_1d(wave_function,scratch_vector)
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
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL cs_diag
      IF(log_main(9)) then
         title='Diagonal Cosine Propagator'
         CALL plot_propagator(cos_diag)
         title='Diagonal Sine Propagator'
         CALL plot_propagator(sin_diag)
      END IF
!
!     Scale wave_function using scratch_vector as scratch
!
!-----------------------------------------------------------------------
      call diagonal_mul(wave_function,scratch_vector)
!-----------------------------------------------------------------------
     IF(log_main(9)) then
        title='First Diagonally Scaled Vector'
        CALL print_psi(wave_function)
     END IF
!
!    Now do the real work, the operation of the off-diagonal propagator
!    on the vector.  
!
     prop_point = 1
     CALL exp_off_diag(wave_function,scratch_vector)
!-----------------------------------------------------------------------
     IF(log_main(9)) then
        title='Off Diagonally Propagated Vector'
        CALL print_psi(wave_function)
     END IF
!
!    Now finish off with the diagonal scale
!-----------------------------------------------------------------------
     call diagonal_mul(wave_function,scratch_vector)
!-----------------------------------------------------------------------
     IF(log_main(9)) then
        title='Second Diagonally Scaled Vector'
        CALL print_psi(wave_function)
     END IF
  ELSE IF (prop_order == 4) THEN
!
!     First Diagonal Scaling by Potential
!     The time is,
!
      tau_loc = p_fac * deltat * .5d0 / hbar
!
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!     First Off-Diagonal Scaling
!     We use the first propagtor
! 
      prop_point = 1
      CALL exp_off_diag(wave_function,scratch_vector)
!
!     Second Diagonal Scaling
!     The time is,
!
      tau_loc = p_fac * deltat / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!     Second Off-Diagonal Scaling
!     No change in which propagator is used.
!
      CALL exp_off_diag(wave_function,scratch_vector)
!
!     Third Diagonal Scaling
!     The time is,
!
      tau_loc = ( 1.d0 - 3.d0 * p_fac ) * deltat * .5d0 / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!     Third Off_diagonal Scaling
!     We need to use the second propagator.
!
      prop_point = 2
      CALL exp_off_diag(wave_function,scratch_vector)
!
!     Fourth Diagonal Scaling
!     Time is same as the third.
!
      call diagonal_mul(wave_function,scratch_vector)
!
!     Fourth Off-Diagonal Scaling
!     BAck to using the first propagator.
!
      prop_point = 1
      CALL exp_off_diag(wave_function,scratch_vector)
!
!     Fifth Diagonal Scaling
!     The time is,
!
      tau_loc = p_fac * deltat / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!     Fifth Off-Diagonal Scaling
!     Still using the first propagator.
!
      CALL exp_off_diag(wave_function,scratch_vector)
!
!     Sixth Diagonal Scaling
!     The time is,
!
      tau_loc = p_fac * deltat * .5d0 / hbar
      CALL cs_diag
      call diagonal_mul(wave_function,scratch_vector)
!
!     OK.  We are done.
! 
  END IF
END SUBROUTINE rsp_1d
