!
MODULE dvr_shared
!***begin prologue     dvr_shared
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           time, finite element dvr, orthogonal polynomial
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            global shared variables for dvr library
!***description        this routine defines the global variables
!***                   and data needed for the dvrlib
!
!***references

!***routines called    
!***end prologue       dvr_shared
  USE input_output
  USE grid_global
  USE potential
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
  CHARACTER (LEN=80), DIMENSION(2)        :: vtyp
  CHARACTER (LEN=80), DIMENSION(30)       :: coord
  CHARACTER (LEN=24)                      :: system
  CHARACTER (LEN=128)                     :: filbec
  CHARACTER (LEN=8)                       :: typke
  INTEGER, DIMENSION(30)                  :: nphy, nglobal, row, nonz
  LOGICAL                                 :: diag, nlse, genpts
!---------------------------------------------------------------------
! Allocated variables for the dvr quantities
!---------------------------------------------------------------------
  REAL*8, DIMENSION(30)                   :: pt_0, pt_n
  TYPE dvr_grid
    REAL*8, DIMENSION(:),   POINTER       :: pt, wt, eigv_0, eigv, v, ang_pot
    REAL*8, DIMENSION(:,:), POINTER       :: f, df, ddf, ke, eigvec_0
    REAL*8, DIMENSION(:,:), POINTER       :: eigvec, p_mom, h
    REAL*8, DIMENSION(:),   POINTER       :: srf_prm
    REAL*8, DIMENSION(:,:), POINTER       :: srf_0, srf
  END TYPE dvr_grid
  TYPE (dvr_grid), DIMENSION(:), ALLOCATABLE &
                                          :: grid
  TYPE buffer
    REAL*8,     DIMENSION(:),   POINTER   :: d, hbuf
    REAL*8,     DIMENSION(:,:), POINTER   :: hibuf
  END TYPE buffer
  TYPE(buffer), DIMENSION(:),                &
                ALLOCATABLE               :: buf    
  REAL*8, DIMENSION(:), ALLOCATABLE       :: v_pert
!
!
END MODULE dvr_shared
