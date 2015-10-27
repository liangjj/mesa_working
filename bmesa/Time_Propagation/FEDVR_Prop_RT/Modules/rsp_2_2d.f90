!*rsp_2_2d.f
!***begin prologue     rsp_2_2d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose            calculate auto correlation function 
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       rsp_2_2d
  SUBROUTINE rsp_2_2d(wave_function,scratch_vector)
     USE dvrprop_global_rt
     USE dvr_shared
     USE exp_off_diagonal
     USE plot_wavefunction
     USE exp_off_diagonal
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),2)             :: wave_function, scratch_vector
!
! Since the potential is, in general, time-dependent, the diagonal
! propagator will be reconstructed at each time-step.
!
!
  IF(log_main(9)) then
     title='Initial vector'
     call print_psi(wave_function)
  END IF
  cos_diag =   cos(v_tot*deltat*.5d0/hbar)  
  sin_diag = - sin(v_tot*deltat*.5d0/hbar)  
  IF(log_main(9)) then
     title='Diagonal Cosine Propagator'
     CALL plot_propagator(cos_diag)
     title='Diagonal Sine Propagator'
     CALL plot_propagator(sin_diag)
  END IF
!
! Scale wave_function using scratch_vector as scratch
!
!-----------------------------------------------------------------------
  scratch_vector = wave_function
  call diagonal_mul(wave_function,scratch_vector,                      &
                    sin_diag,cos_diag,n3d)
!-----------------------------------------------------------------------
  IF(log_main(9)) then
     title='First Diagonally Scaled Vector'
     CALL print_psi(wave_function)
  END IF
!
! Now do the real work, the operation of the off-diagonal propagator
! on the vector.  
!----------------------------------------------------------------------- 
  CALL exp_off_diag(wave_function,scratch_vector)
!-----------------------------------------------------------------------
  IF(log_main(9)) then
     title='Off Diagonally Propagated Vector'
     CALL print_psi(wave_function)
  END IF
!
! Now finish off with the diagonal scale
!-----------------------------------------------------------------------
  scratch_vector = wave_function
  call diagonal_mul(wave_function,scratch_vector,                      &
                    sin_diag,cos_diag,n3d)
!-----------------------------------------------------------------------
  IF(log_main(9)) then
     title='Second Diagonally Scaled Vector'
     CALL print_psi(wave_function)
  END IF
END SUBROUTINE rsp_2_2d
