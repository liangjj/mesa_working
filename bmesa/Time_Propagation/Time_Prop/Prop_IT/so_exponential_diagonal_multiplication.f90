
!***********************************************************************
                           MODULE so_exponential_diagonal_multiplication
                           INTERFACE so_exp_diagonal_multiply
                    MODULE PROCEDURE so_exp_diagonal_mul_d,              &
                                     so_exp_diagonal_mul_z
                       END INTERFACE so_exp_diagonal_multiply
!***********************************************************************
                           CONTAINS
!***********************************************************************
!***********************************************************************
!***********************************************************************
!deck so_exp_diagonal_mul_d
!***begin prologue     so_exp_diagonal_mul_d
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
!***end prologue       so_exp_diagonal_mul_d
!
  SUBROUTINE so_exp_diagonal_mul_d(vector_in_out)
  USE dvrprop_global
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)                         :: vector_in_out
  vector_in_out(:) = exp_diag_d(:) * vector_in_out(:)  
END SUBROUTINE so_exp_diagonal_mul_d
!***********************************************************************
!***********************************************************************
!deck so_exp_diagonal_mul_z
!***begin prologue     exp_diagonal_mul_z
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
!***end prologue       so_exp_diagonal_mul_z
!
  SUBROUTINE so_exp_diagonal_mul_z(vector_in_out)
  USE dvrprop_global
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d)                       :: vector_in_out
  vector_in_out(:) = exp_diag_z(:) * vector_in_out(:)  
END SUBROUTINE so_exp_diagonal_mul_z
!***********************************************************************
!***********************************************************************
END MODULE so_exponential_diagonal_multiplication
