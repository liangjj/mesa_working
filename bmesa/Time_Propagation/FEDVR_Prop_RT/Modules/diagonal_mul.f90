!deck diagonal_mul
!***begin prologue     diagonal_mul
!***date written       040707   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            diagonal scaling of propagated vector
!***
!***description        
!                      
!                      
!***references
!***routines called    
!***end prologue       diagonal_mul
!
  SUBROUTINE diagonal_mul(vector_in_out,vector_temp)
  USE dvrprop_global_rt
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d,2)                       :: vector_in_out
  REAL*8, DIMENSION(n3d,2)                       :: vector_temp
  vector_temp(:,1) = cos_diag(:) * vector_in_out(:,1)       &
                                 -                          &
                     sin_diag(:) * vector_in_out(:,2)  
  vector_temp(:,2) = sin_diag(:) * vector_in_out(:,1)       &
                                 +                          &
                     cos_diag(:) * vector_in_out(:,2)  
  vector_in_out = vector_temp
END SUBROUTINE diagonal_mul
