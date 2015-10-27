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
!                   Keywords for Print Options
!
  CHARACTER (LEN=80), DIMENSION(7)       :: pr_prp
  DATA pr_prp  / 'trial','h-on-vectors','small-matrix','eigenvalues', &
                 'eigenvectors','solution','none' /
!
!
END MODULE lanczos_prnt
