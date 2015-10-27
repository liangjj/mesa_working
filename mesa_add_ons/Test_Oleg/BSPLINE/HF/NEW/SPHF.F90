!======================================================================
  PROGRAM Spline_HF
!======================================================================
!
!  This program computes the radial functions for simple Hartree-Fock
!  cases.
!
!   SUBROUTINE called:
!       get_case
!       define_grid
!       define_spline
!       define_atomic
!       obtain_estimates
!       optimize_orbitals
!       properties
!
!   Date: 01/21/1999
!                                                                
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    REAL(KIND=8) :: z

    ! .. get data about the problem to be solved
    CALL get_case(z)

    ! .. sets up positions of grid points: t
    CALL define_grid (z)

    ! .. initializes the values of the spline and its derivatives
    ! .. and evaluates the spline arrays (operators in spline basis)
    ! .. which are defined in the MODULES spline_grid and spline_galerkin
    CALL define_spline(z)

    ! ... initialize that atomic physics environment
    CALL define_atomic

    ! .. obtain initial estimates for orbitals
    CALL obtain_estimates

    ! .. optimize the orbitals for a stationary solution
    CALL optimize_orbitals

    ! .. compute atomic properties
    CALL properties

  END PROGRAM slater
