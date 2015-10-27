!
MODULE dvr_prnt
!***begin prologue     dvr_prnt
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           time, finite element dvr, orthogonal polynomial
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            global print variables for dvr library
!***description        this routine defines the global print
!***                   variables
!
!***references

!***routines called    
!***end prologue       dvr_prnt
  USE input_output

  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                   Keywords for Print Options
!
   CHARACTER (LEN=80), DIMENSION(12) :: prnkey
   CHARACTER (LEN=80)                :: secprn
   DATA prnkey /'sector-points','sector-polynomials', &
                'sector-matrix','global-points','global-polynomials', &
                'potential','global-matrix-elements','hamiltonian',  &
                'eigenvalues','eigenvectors','amplitudes', &
                'none' /
   DATA secprn / 'none' /
!
!
  CHARACTER (LEN=80),DIMENSION(12)      :: prloc
  LOGICAL, DIMENSION(12)                :: prn
!
!
END MODULE dvr_prnt
