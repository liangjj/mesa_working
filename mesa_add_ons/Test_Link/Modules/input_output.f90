!
MODULE input_output
!***begin prologue     input_output
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source             
!***purpose            
!
!***references

!***routines called    
!***end prologue       input_output

  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                   Some Constants and Units
!
   INTEGER                                       :: inp=8
   INTEGER                                       :: iout=9
   INTEGER                                       :: rows_to_print
   INTEGER                                       :: columns_to_print
   INTEGER                                       :: eigenvectors_to_print
   LOGICAL                                       :: print_parameter
   CHARACTER(LEN=8), DIMENSION(:), ALLOCATABLE   :: rowlab, collab

!
!
END MODULE input_output
