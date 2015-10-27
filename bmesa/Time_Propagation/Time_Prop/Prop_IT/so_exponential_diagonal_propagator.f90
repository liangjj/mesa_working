!***********************************************************************
                           MODULE so_exponential_diagonal_propagator
                         INTERFACE so_exp_diagonal_propagator
                  MODULE PROCEDURE so_exp_diag_prop_d,                   &
                                   so_exp_diag_prop_z
                         END INTERFACE so_exp_diagonal_propagator
!***********************************************************************
                           CONTAINS
!***********************************************************************
!***********************************************************************
!***********************************************************************
!deck so_exp_diag_prop_d
!***begin prologue     so_exp_diag_prop_d
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
!***end prologue       so_exp_diag_prop_d
  SUBROUTINE so_exp_diag_prop_d(exp_diag)
  USE dvrprop_global
  USE dvr_shared
  USE arnoldi_global
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)                  :: exp_diag
  exp_diag =   exp(- v_tot * tau_loc)  
END SUBROUTINE so_exp_diag_prop_d
!***********************************************************************
!***********************************************************************
!deck so_exp_diag_prop_z
!***begin prologue     so_exp_diag_prop_z
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
!***end prologue       so_exp_diag_prop_z
  SUBROUTINE so_exp_diag_prop_z(exp_diag)
  USE dvrprop_global
  USE dvr_shared
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d)                  :: exp_diag
  exp_diag =   exp( - eye * v_tot * tau_loc)  
END SUBROUTINE so_exp_diag_prop_z
!***********************************************************************
!***********************************************************************
                     END MODULE so_exponential_diagonal_propagator
