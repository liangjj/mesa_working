!
MODULE dvr_global
!***begin prologue     dvr_global
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           time, finite element dvr, orthogonal polynomial
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            global variables for dvr library
!***description        this routine defines the global variables
!***                   and data needed for the dvrlib
!
!***references

!***routines called    
!***end prologue       dvr_global
  USE input_output
  USE grid_global
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
  CHARACTER (LEN=4)                     :: parity
  INTEGER, PARAMETER                    :: maxreg=5000
  INTEGER, DIMENSION(maxreg)            :: n
  INTEGER, DIMENSION(maxreg)            :: npt
  INTEGER, DIMENSION(maxreg)            :: nrq
  INTEGER                               :: spdim
  INTEGER                               :: nfix, bcl, bcr
  INTEGER                               :: nreg, l_orb_max, l_max 
  INTEGER                               :: l_val, m_val 
  INTEGER                               :: nphy_tri 
  REAL*8,  DIMENSION(maxreg+1)          :: edge
  REAL*8                                :: box, deltax
  LOGICAL, DIMENSION(2)                 :: drop, fix
  LOGICAL                               :: skip
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: nfun_reg
  INTEGER, DIMENSION(:), ALLOCATABLE    :: num_reg
  REAL*4,  DIMENSION(20)                :: del, tim
!
!
END MODULE dvr_global
