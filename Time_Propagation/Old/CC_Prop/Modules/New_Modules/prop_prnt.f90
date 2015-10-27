!
MODULE prop_prnt
!***begin prologue     prop_prnt
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
!***end prologue       prop_prnt
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                   Keywords for Print Options
!
  CHARACTER (LEN=80), DIMENSION(9)       :: pr_prp
  CHARACTER (LEN=80), DIMENSION(9)       :: pr_main
  DATA pr_prp  / 'trials','orthogonalize-trials',  &
                 'h-on-vectors','small-matrix','eigenvalues', &
                 'eigenvectors','overlaps','solution','none' /
  DATA pr_main / 'eigenvalues','eigenvectors',  &
                 'potential','initial-state', &
                 'non-linear-potential','h-on-initial-state', &
                 'trial-vectors','solution','none' /
!
!
END MODULE prop_prnt
