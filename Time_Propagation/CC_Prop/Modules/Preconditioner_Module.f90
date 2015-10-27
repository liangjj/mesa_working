!***********************************************************************
                           MODULE Preconditioner_Module
                           USE dvrprop_global
                           USE Iterative_Global
                           USE Pack_Global
                           USE Pack_Hamiltonian_Module
                           USE Matrix_Utility_Module
                           IMPLICIT NONE
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
                          INTERFACE Pack_Triangle
                   MODULE PROCEDURE Pack_Triangle_d,                          &
                                    Pack_Triangle_z
                          END INTERFACE Pack_Triangle
!
                          INTERFACE Triangular_Solve
                   MODULE PROCEDURE Triangular_Solve_d,                       &
                                    Triangular_Solve_z
                          END INTERFACE Triangular_Solve
!
                          INTERFACE Packed_Solve
                   MODULE PROCEDURE Packed_Solve_d,                           &
                                    Packed_Solve_z
                          END INTERFACE Packed_Solve
!
                          INTERFACE Upper_to_Lower
                   MODULE PROCEDURE Upper_to_Lower_d,                         &
                                    Upper_to_Lower_z
                          END INTERFACE Upper_to_Lower
!
                          INTERFACE Full_to_Upper
                   MODULE PROCEDURE Full_to_Upper_d,                          &
                                    Full_to_Upper_z
                          END INTERFACE Full_to_Upper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           Contains
!***********************************************************************
!***********************************************************************
!deck Full_to_Upper_d
!***begin prologue     Full_to_Upper_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             cholesky decomposition of real, symmetric, positive
!***                   definite matrix
!***purpose            
!***
!***references
!***routines called
!***end prologue       Full_to_Upper_d
!
  SUBROUTINE Full_to_Upper_d (matrix,upper)
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)              :: matrix
  REAL*8, DIMENSION(:)                :: upper
  INTEGER                             :: i
  INTEGER                             :: j
  INTEGER                             :: count
!
  count = 0
  DO i=1,n3d
     DO j=1, i
        count = count + 1
        upper(count) = matrix(i,j)
     END DO
  END DO
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Full_to_Upper_d
!***********************************************************************
!***********************************************************************
!deck Full_to_Upper_z
!***begin prologue     Full_to_Upper_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             cholesky decomposition of real, symmetric, positive
!***                   definite matrix
!***purpose            
!***
!***references
!***routines called
!***end prologue       Full_to_Upper_z
!
  SUBROUTINE Full_to_Upper_z (matrix,upper)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:,:)          :: matrix
  COMPLEX*16, DIMENSION(:)            :: upper
  INTEGER                             :: i
  INTEGER                             :: j
  INTEGER                             :: count
!
  count = 0
  DO i=1,n3d
     DO j=1, i
        count = count + 1
        upper(count) = matrix(i,j)
     END DO
  END DO
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Full_to_Upper_z
!***********************************************************************
!***********************************************************************
!deck Pack_Triangle_d
!***begin prologue     Pack_Triangle_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             cholesky decomposition of real, symmetric, positive
!***                   definite matrix
!***purpose            
!***
!***references
!***routines called
!***end prologue       Pack_Triangle_d
!
  SUBROUTINE Pack_Triangle_d (type,triangle,packed_columns,non_zero_columns, &
                              row_index,number,matrix_diagonal)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                :: triangle
  REAL*8, DIMENSION(:)                :: matrix_diagonal
  REAL*8, DIMENSION(:)                :: packed_columns
  INTEGER, DIMENSION(:)               :: non_zero_columns
  INTEGER, DIMENSION(:)               :: row_index
  INTEGER                             :: number
  CHARACTER(LEN=*)                    :: type 
!
!
! Pack a triangle Matrix to efficiently solve the linear equations
!
  IF ( type /= 'lower') THEN
       CALL h_pack_upper(triangle,                                           &
                         matrix_diagonal,                                    &
                         packed_columns,                                     &
                         non_zero_columns,                                   &
                         row_index,                                          &
                         number,n3d,type)
       Call Read_and_Write_Column_Packed_Matrices (packed_columns,           &
                                                   non_zero_columns,         &
                                                   row_index,                &
                                                  'write',                   &
                                                   type,                     &
                                                   number,                   &
                                                   matrix_diagonal )
  ELSE    
      CALL h_pack_lower(triangle,                                            &
                        matrix_diagonal,                                     &
                        packed_columns,                                      &
                        non_zero_columns,                                    &
                        row_index,                                           &
                        number,n3d,type)
      Call Read_and_Write_Column_Packed_Matrices (packed_columns,            &
                                                  non_zero_columns,          &
                                                  row_index,                 &
                                                 'write',                    &
                                                  type,                      &
                                                  number )
  END IF
1 Format(/,5x,'Column = ',i4)
2  Format( (15x,5f10.5) )
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Pack_Triangle_d
!***********************************************************************
!***********************************************************************
!***begin prologue     Pack_Triangle_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source             
!***                   
!***purpose            
!***
!***references
!***routines called
!***end prologue       Pack_Triangle_z
!
  SUBROUTINE Pack_Triangle_z (type,triangle,packed_columns,non_zero_columns, &
                              row_index,number,matrix_diagonal)
  COMPLEX*16, DIMENSION(:)            :: triangle
  COMPLEX*16, DIMENSION(:), OPTIONAL  :: matrix_diagonal
  COMPLEX*16, DIMENSION(:)            :: packed_columns
  INTEGER, DIMENSION(:)               :: non_zero_columns
  INTEGER, DIMENSION(:)               :: row_index
  INTEGER                             :: number
  CHARACTER(LEN=*)                    :: type 
!
! Pack the Cholesky factors in upper and lower form to efficiently
! solve the linear equations
! 
  IF ( type /= 'lower') THEN
       CALL h_pack_upper(triangle,                                           &
                         matrix_diagonal,                                    &
                         packed_columns,                                     &
                         non_zero_columns,                                   &
                         row_index,                                          &
                         number,n3d,type)
       Call Read_and_Write_Column_Packed_Matrices (packed_columns,           &
                                                   non_zero_columns,         &
                                                   row_index,                &
                                                  'write',                   &
                                                   type,                     &
                                                   number,                   &
                                                   matrix_diagonal )
  ELSE    
      CALL h_pack_lower(triangle,                                            &
                        matrix_diagonal,                                     &
                        packed_columns,                                      &
                        non_zero_columns,                                    &
                        row_index,                                           &
                        number,n3d,type)
      Call Read_and_Write_Column_Packed_Matrices (packed_columns,            &
                                                  non_zero_columns,          &
                                                  row_index,                 &
                                                 'write',                    &
                                                  type,                      &
                                                  number )
  END IF
!
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Pack_Triangle_z
!***********************************************************************
!**********************************************************************
!deck Triangular_Solve_d
!***begin prologue     Triangular_Solve_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Triangular_Solve_d
!
  SUBROUTINE Triangular_Solve_d (solution,right_hand_side,trans,it)
  IMPLICIT NONE
  REAL*8,  DIMENSION(:)                 :: solution
  REAL*8,  DIMENSION(:)                 :: right_hand_side
  CHARACTER(LEN=1)                      :: trans
  CHARACTER(LEN=3)                      :: itoc
  INTEGER                               :: it  
!
  solution = right_hand_side
  Call dtpsv('u',trans,'n',n3d,upper_d,solution,1)
  IF(print_internal_matrices) THEN
     title='Solution Vector Iteration = '//itoc(it)//' t = '//trans
     call prntfmn(title,solution,n3d,1,n3d,1,iout,'e')
  END IF
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Triangular_Solve_d
!***********************************************************************
!***********************************************************************
!deck Triangular_Solve_z
!***begin prologue     Triangular_Solve_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Triangular_Solve_z
!
  SUBROUTINE Triangular_Solve_z (solution,right_hand_side,trans,it)
  IMPLICIT NONE
  COMPLEX*16,  DIMENSION(:)             :: solution
  COMPLEX*16,  DIMENSION(:)             :: right_hand_side
  CHARACTER(LEN=1)                      :: trans
  CHARACTER(LEN=3)                      :: itoc
  INTEGER                               :: it  
!
  solution = right_hand_side
  Call ztpsv('u',trans,'n',n3d,upper_z,solution,1)
  IF(print_internal_matrices) THEN
     title='Solution Vector Iteration = '//itoc(it)//' t = '//trans
     call prntcmn(title,solution,n3d,1,n3d,1,iout,'e')
  END IF
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Triangular_Solve_z
!***********************************************************************
!***********************************************************************
!deck Packed_Solve_d
!***begin prologue     Packed_Solve_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Packed_Solve_d
!
  SUBROUTINE Packed_Solve_d (solution,right_hand_side,non_zero_columns,row_index,         &
                             packed_columns,matrix_diagonal,number,trans,it)
  IMPLICIT NONE
  REAL*8,  DIMENSION(:)                 :: solution
  REAL*8,  DIMENSION(:)                 :: right_hand_side
  REAL*8,  DIMENSION(:)                 :: packed_columns
  REAL*8,  DIMENSION(:)                 :: matrix_diagonal
  INTEGER, DIMENSION(:)                 :: non_zero_columns
  INTEGER, DIMENSION(:)                 :: row_index
  INTEGER                               :: number
  CHARACTER(LEN=1)                      :: trans
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: ij
  INTEGER                               :: it  
  CHARACTER(LEN=3)                      :: itoc
  INTEGER                               :: count
!
  solution = right_hand_side
  IF (trans == 'n') THEN
!
!     This is going to be a solve in the backward direction using
!     the packed upper triangle.
!
      count = number
      DO i=n3d, 1, -1
!
!        Step one is to solve for the diagonal element
!
         solution(i) = solution(i) / matrix_diagonal(i)  
!
!        Now we must eliminate this unknown from all of the other equations.
!  
         DO j=1,non_zero_columns(i)
            ij=row_index(count)
            solution(ij) = solution(ij)                                  &
                                  -                                      &
                           packed_columns(count) * solution(i)
            count = count - 1
         END DO
      END DO
  ELSE IF (trans == 't') THEN
!
!      This is going to be a solve in the forward direction using
!      the lower triangle.
!
         count = 1
         DO i=1, n3d
!
!           Step one is to solve for the diagonal element
! 
            solution(i) = solution(i) / matrix_diagonal(i)  
!
!           Now we must eliminate this unknown from all of the other equations.
!  
            DO j=1,non_zero_columns(i)
               ij=row_index(count)
               solution(ij) = solution(ij)                                  &
                                     -                                      &
                              packed_columns(count) * solution(i)
               count = count + 1
            END DO
         END DO
  END IF
  IF(print_internal_matrices) THEN
     title='Solution Vector Iteration = '//itoc(it)//' t = '//trans
     call prntfmn(title,solution,n3d,1,n3d,1,iout,'e')
  END IF
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Packed_Solve_d
!***********************************************************************
!***********************************************************************
!deck Packed_Solve_z
!***begin prologue     Packed_Solve_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Packed_Solve_z
!
  SUBROUTINE Packed_Solve_z (solution,right_hand_side,non_zero_columns,row_index,         &
                             packed_columns,matrix_diagonal,number,trans,it)
  IMPLICIT NONE
  COMPLEX*16,  DIMENSION(:)             :: solution
  COMPLEX*16,  DIMENSION(:)             :: right_hand_side
  COMPLEX*16,  DIMENSION(:)             :: packed_columns
  COMPLEX*16,  DIMENSION(:)             :: matrix_diagonal
  COMPLEX*16                            :: conjg
  INTEGER, DIMENSION(:)                 :: non_zero_columns
  INTEGER, DIMENSION(:)                 :: row_index
  INTEGER                               :: number
  CHARACTER(LEN=1)                      :: trans
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: ij
  CHARACTER(LEN=3)                      :: itoc
  INTEGER                               :: it  
  INTEGER                               :: count
  INTEGER                               :: ipnt
!
  solution = right_hand_side
  IF (trans == 'n') THEN
!
!     This is going to be a solve in the backward direction
!
      count = number
      DO i=n3d, 1, -1
!        Step one is to solve for the diagonal element
!
         solution(i) = solution(i) / matrix_diagonal(i)  
!
!         Now we must eliminate this unknown from all of the other equations.
!  
         DO j=1,non_zero_columns(i)
            ij=row_index(count)
            solution(ij) = solution(ij)                                  &
                                  -                                      &
                           packed_columns(count) * solution(i)
            count = count - 1
         END DO
      END DO
  ELSE IF (trans == 'c') THEN
!
!     This is going to be a solve in the forward direction
!
      count = 1
      DO i=1, n3d
!
!        Step one is to solve for the diagonal element
!
         solution(i) = solution(i) / matrix_diagonal(i)  
!        Now we must eliminate this unknown from all of the other equations.
!  
         DO j=1,non_zero_columns(i)
            ij=row_index(count)
            solution(ij) = solution(ij)                                  &
                                  -                                      &
                           conjg(packed_columns(count)) * solution(i)
            count = count + 1
         END DO
      END DO
  END IF
  IF(print_internal_matrices) THEN
     title='Solution Vector Iteration = '//itoc(it)//' t = '//trans
     call prntcmn(title,solution,n3d,1,n3d,1,iout,'e')
  END IF
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Packed_Solve_z
!***********************************************************************
!***********************************************************************
!deck Upper_to_Lower_d
!***begin prologue     Upper_to_Lower_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Upper_to_Lower_d
!
  SUBROUTINE Upper_to_Lower_d (upper,lower)
  IMPLICIT NONE
  REAL*8,  DIMENSION(:)                 :: upper
  REAL*8,  DIMENSION(:)                 :: lower
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: ij
  INTEGER                               :: count
!
!
  count = 0
  DO j = 1,n3d
     DO i = j, n3d
        ij = i*(i-1)/2 + j
        count = count + 1
        lower(count) = upper(ij)
     END DO
  END DO
!
!
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Upper_to_Lower_d
!***********************************************************************
!***********************************************************************
!deck Upper_to_Lower_z
!***begin prologue     Upper_to_Lower_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Upper_to_Lower_z
!
  SUBROUTINE Upper_to_Lower_z (upper,lower)
  IMPLICIT NONE
  COMPLEX*16,  DIMENSION(:)             :: upper
  COMPLEX*16,  DIMENSION(:)             :: lower
  COMPLEX*16                            :: conjg
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: ij
  INTEGER                               :: count
!
!
  count = 0
  DO j =1,n3d
     DO i=j,n3d
        ij = i*(i-1)/2 + j
        count = count + 1
        lower(count) = conjg(upper(ij))
     END DO
  END DO
!
!
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Upper_to_Lower_z
!***********************************************************************
!***********************************************************************
  END  MODULE Preconditioner_Module
!***********************************************************************
!***********************************************************************



