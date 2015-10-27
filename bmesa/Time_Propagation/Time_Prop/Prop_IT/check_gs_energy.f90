!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck check_gs_energy.f
!***begin prologue     check_gs_energy
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check energy one dimension wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_gs_energy
  SUBROUTINE check_gs_energy(v,hv,e)
  USE dvrprop_global
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)             :: v
  REAL*8, DIMENSION(n3d)             :: hv
  REAL*8                             :: ddot
  REAL*8                             :: e
!
  e = ddot(n3d,v,1,hv,1)
END SUBROUTINE check_gs_energy
!***********************************************************************
