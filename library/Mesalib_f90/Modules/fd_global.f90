!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE fd_global
!***begin prologue     fd_global
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           finite difference
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            global variables for finite difference routines
!***references

!***routines called    
!***end prologue       fd_global
  USE input_output
  USE grid_global
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
  INTEGER                               :: nstep, ndiff
  REAL*8                                :: del, dscale
  REAL*8, DIMENSION(4)                  :: d
  REAL*8, DIMENSION(2)                  :: edge
!
!
END MODULE fd_global
