!deck real_space_propagator_2_order
!***begin prologue     real_space_propagator_2_order
!***date written       031122   (yymmdd)
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
!***end prologue       real_space_propagator_2_order
!
  SUBROUTINE real_space_propagator_2_order
  USE dvrprop_global
  USE dvr_shared
  USE exp_off_diagonal
  IMPLICIT NONE
  INTEGER                             :: i
!
! Since the potential is, in general, time-dependent, the diagonal
! propagator will be reconstructed at each time-step.
!
!
  IF(log_main(9)) then
     title='Initial vector'
     call plot_wavefunction(psi_08)
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
! Scale psi_08 using v_scr_08 as scratch
!
!-----------------------------------------------------------------------
  v_scr_08 = psi_08
  psi_08(:,1) = cos_diag(:) * v_scr_08(:,1) - sin_diag(:) * v_scr_08(:,2)   
  psi_08(:,2) = sin_diag(:) * v_scr_08(:,1) + cos_diag(:) * v_scr_08(:,2)  
!-----------------------------------------------------------------------
  IF(log_main(9)) then
     title='First Diagonally Scaled Vector'
     CALL plot_wavefunction(psi_08)
  END IF
!
! Now do the real work, the operation of the off-diagonal propagator
! on the vector.  
!----------------------------------------------------------------------- 
  CALL exp_off_diag
!-----------------------------------------------------------------------
  IF(log_main(9)) then
     title='Off Diagonally Propagated Vector'
     CALL plot_wavefunction(psi_08)
  END IF
!
! Now finish off with the diagonal scale
!-----------------------------------------------------------------------
  v_scr_08 = psi_08
  psi_08(:,1) = cos_diag(:) * v_scr_08(:,1) - sin_diag(:) * v_scr_08(:,2)   
  psi_08(:,2) = sin_diag(:) * v_scr_08(:,1) + cos_diag(:) * v_scr_08(:,2)  
!-----------------------------------------------------------------------
  IF(log_main(9)) then
     title='Second Diagonally Scaled Vector'
     CALL plot_wavefunction(psi_08)
  END IF
END SUBROUTINE real_space_propagator_2_order
