!
MODULE lanczos_global
!***begin prologue     lanczos_global
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lanczos, diagonalization, hermitian, eigenvalues
!***author             schneider, b. i.(nsf)
!***source             lancz
!***purpose            global variables for lanczos diagonalization
!***description        this routine defines the global variables
!***                   and data needed for the dvrlib
!
!***references

!***routines called    
!***end prologue       lanczos_global
  USE io
  USE prop_global
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
  REAL*8,            DIMENSION(:),                                &
                     ALLOCATABLE      :: a_lanczos, b_lanczos
  REAL*8,            DIMENSION(:,:),                              &
                     ALLOCATABLE      :: b
  REAL*8,            DIMENSION(:),                                &
                     ALLOCATABLE      :: bwrk
  COMPLEX*16,        DIMENSION(:,:),                              &
                     ALLOCATABLE      :: vec, u
  COMPLEX*16,        DIMENSION(:),                                &
                     ALLOCATABLE      :: hvec
  REAL*4,            DIMENSION(10)    :: eltim, del_t 
  COMPLEX*16,        DIMENSION(:),                                &
                     ALLOCATABLE      :: psi0, chi, soln_0      
  COMPLEX*16,        DIMENSION(:),                                &
                     ALLOCATABLE      :: work, vscr
  REAL*8,            DIMENSION(:),                                &
                     ALLOCATABLE      :: eig
!
!
END MODULE lanczos_global
