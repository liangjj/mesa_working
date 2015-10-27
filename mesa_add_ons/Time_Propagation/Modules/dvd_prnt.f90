!
MODULE dvd_prnt
!***begin prologue     dvd_prnt
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           davidson, eigenvalue
!***author             schneider, b. i.(nsf)
!***source             dvdlib
!***purpose            global print variables for davidson library
!***description        this routine defines the global print
!***                   variables
!
!***references

!***routines called    
!***end prologue       dvd_prnt
  USE input_output

  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                   Keywords for Print Options
!
   CHARACTER (LEN=80), DIMENSION(12) :: prn_dvd
   CHARACTER (LEN=80), DIMENSION(12) :: prn_dvd_loc
   LOGICAL, DIMENSION(12)            :: log_dvd
   DATA prn_dvd / 'trials','orthogonalized-trials',               &
                  'h-on-vectors','small-matrix','eigenvalues',    &
                  'iteration-information','transformed-vectors',  &
                  'transformed-h-on-vectors','residuals',         &
                  'new-raw-vectors','overlaps','all' /          

!
!
END MODULE dvd_prnt
