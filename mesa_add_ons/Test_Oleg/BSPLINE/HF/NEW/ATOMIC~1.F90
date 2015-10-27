!=======================================================================
  MODULE atomic_state 
!=======================================================================
!   This module defines the parameters for the problem to be solved
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    SAVE

    ! atom variables
    REAL(KIND=8) ::z
    INTEGER :: nclosd, nwf, nit
    CHARACTER(LEN=2) :: atom
    CHARACTER(LEN=3) :: term
    CHARACTER(LEN=32) :: config

    ! orbital variables
    CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: el
    INTEGER, DIMENSION(:), ALLOCATABLE :: n, l, max, meth, ind
    INTEGER :: lmax
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: p
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: az, dpm, s, sum

    ! energy expression variables
    INTEGER :: kmax
    INTEGER, DIMENSION(:), ALLOCATABLE :: ijptr
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: coef
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: e

    CONTAINS

    !===================================================================
      SUBROUTINE allocate_atomic_state
    !===================================================================
    !   This program allocates arrays associated with number of orbitals
    !-------------------------------------------------------------------
        IMPLICIT NONE
        
	ALLOCATE( n(nwf), l(nwf), max(nwf), ind(nwf), s(nwf), meth(nwf))
	ALLOCATE( az(nwf), dpm(nwf), s(nwf), sum(nwf) )
        ALLOCATE( e(nwf,nwf) )

      END SUBROUTINE allocate_atomic_state

    !===================================================================
      SUBROUTINE allocate_orbital_array
    !===================================================================
    !   This program allocates arrays associated with number of orbitals
    !   and the spline expansion
    !-------------------------------------------------------------------
        USE spline_param
        IMPLICIT NONE
        
	ALLOCATE( p(ns,nwf) )
        ALLOCATE( ijptr(nwf-nclosd, nwf-nclosd) )
     
      END SUBROUTINE allocate_orbital_array
