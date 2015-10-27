!deck exp_diag_prop.f
!***begin prologue     exp_diag_prop
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            time-dependent wave function calculation
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       exp_diag_prop
  SUBROUTINE exp_diag_prop
     USE dvrprop_global_it
     USE dvr_shared
  IMPLICIT NONE
  exp_diag =   exp(-v_tot * tau_loc)  
END SUBROUTINE exp_diag_prop
