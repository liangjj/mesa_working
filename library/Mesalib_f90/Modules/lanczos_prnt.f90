!
MODULE lanczos_prnt
!***begin prologue     lanczos_prnt
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           propagation
!***author             schneider, b. i.(nsf)
!***source             proplib
!***purpose            global print variables for propagation
!***description        this routine defines the global print
!***                   variables
!
!***references

!***routines called    
!***end prologue       lanczos_prnt
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                   Keywords for Lanczos Print Options
!
  CHARACTER (LEN=80), DIMENSION(10)       :: print_lanczos
  LOGICAL,            DIMENSION(10)       :: log_lanczos
  DATA print_lanczos  / 'recursion_coefficients','lanczos_vectors',           &
                        'h_on_vector','overlap_matrix', 'small_matrices',     &
                        'eigenvalues', 'eigenvectors','convergence_tests',    &
                        'solution', 'none' /
!
!
END MODULE lanczos_prnt
