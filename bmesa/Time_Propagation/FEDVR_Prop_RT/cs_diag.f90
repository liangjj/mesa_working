!deck cs_diag.f
!***begin prologue     cs_diag
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            time-dependent wave function calculation
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       cs_diag
  SUBROUTINE cs_diag
     USE dvrprop_global_rt
     USE dvr_shared
  IMPLICIT NONE
  cos_diag =   cos(v_tot * tau_loc)  
  sin_diag = - sin(v_tot * tau_loc)  
END SUBROUTINE cs_diag
