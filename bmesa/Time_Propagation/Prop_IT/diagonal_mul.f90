!deck diagonal_mul
!***begin prologue     diagonal_mul
!***date written       040707   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            diagonal scaling of propagated vector by exponential.
!***
!***description        
!                      
!                      
!***references
!***routines called    
!***end prologue       diagonal_mul
!
  SUBROUTINE diagonal_mul(vector_in_out)
  USE dvrprop_global_it
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)                         :: vector_in_out
  vector_in_out = exp_diag * vector_in_out 
END SUBROUTINE diagonal_mul
