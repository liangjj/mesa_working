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
  CHARACTER (LEN=24)                      :: kinetic_energy_type
  INTEGER, DIMENSION(30)                  :: nphy, nglobal, row, nonz
  LOGICAL                                 :: diag, nlse, genpts
!---------------------------------------------------------------------
! Allocated variables for the dvr quantities
!---------------------------------------------------------------------
  REAL*8, DIMENSION(30)                       :: pt_0, pt_n
  TYPE dvr_grid
    REAL*8, DIMENSION(:),   ALLOCATABLE       :: pt, wt, eigv_0, eigv, v, ang_pot
    REAL*8, DIMENSION(:,:), ALLOCATABLE       :: f, df, ddf, ke, eigvec_0
    REAL*8, DIMENSION(:,:), ALLOCATABLE       :: eigvec, p_mom, h
    REAL*8, DIMENSION(:),   ALLOCATABLE       :: srf_prm
    REAL*8, DIMENSION(:,:), ALLOCATABLE       :: srf_0, srf
  END TYPE dvr_grid
  TYPE (dvr_grid), DIMENSION(:), ALLOCATABLE &
                                              :: grid
  TYPE pt_wt
       REAL*8, DIMENSION(:), ALLOCATABLE      :: qr
       REAL*8, DIMENSION(:), ALLOCATABLE      :: wtr
       REAL*8, DIMENSION(:), ALLOCATABLE      :: inv_sqrt_wtr
  END TYPE pt_wt

  TYPE poly
       REAL*8, DIMENSION(:,:), ALLOCATABLE    :: pr
       REAL*8, DIMENSION(:,:), ALLOCATABLE    :: dpr
       REAL*8, DIMENSION(:,:), ALLOCATABLE    :: ddpr
  END TYPE poly

  TYPE mat
       REAL*8, DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL*8, DIMENSION(:), ALLOCATABLE      :: vr
       REAL*8, DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE mat

  TYPE regional_grid
       INTEGER                                :: num_reg
       INTEGER, DIMENSION(:), ALLOCATABLE     :: num_pts_reg
       TYPE(pt_wt), DIMENSION(:),                                     &
                       ALLOCATABLE            :: reg_pt_wt
       TYPE(poly), DIMENSION(:,:),                                    &
                       ALLOCATABLE            :: reg_poly
       TYPE(mat), DIMENSION(:,:),                                     &
                       ALLOCATABLE            :: reg_mat
  END TYPE regional_grid

  TYPE(regional_grid), DIMENSION(:),                                  &
                       ALLOCATABLE            :: reg_grid

  TYPE buffer
    REAL*8,     DIMENSION(:), ALLOCATABLE                   &
                                              :: d, hbuf
    REAL*8,     DIMENSION(:,:), ALLOCATABLE                 &
                                              :: hibuf
  END TYPE buffer
  TYPE(buffer), DIMENSION(:),                               &
                ALLOCATABLE                   :: buf    
  REAL*8, DIMENSION(:), ALLOCATABLE           :: v_pert
!
!
END MODULE dvr_shared
