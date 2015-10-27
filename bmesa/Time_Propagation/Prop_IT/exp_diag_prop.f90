!deck exp_diag_prop.f
!***begin prologue     exp_diag_prop
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate exponential time-dependent scaling
!***                   for a local potential, v_tot.
!***description        In the split operator approach, local potentials,
!***                   which may be time dependent are split off before
!***                   treating the kinetic energy operators.  These only
!***                   scale the vectors and do not require the more complex
!***                   even/odd splitting of the kinetic energy.
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       exp_diag_prop
  SUBROUTINE exp_diag_prop
     USE dvrprop_global_it
     USE dvr_shared
     USE arnoldi_global_it
  IMPLICIT NONE
  exp_diag =   exp(- v_tot * tau_loc)  
END SUBROUTINE exp_diag_prop
