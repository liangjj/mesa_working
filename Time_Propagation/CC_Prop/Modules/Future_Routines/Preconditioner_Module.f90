!***********************************************************************
                           MODULE Preconditioner_Module
                           USE dvrprop_global
                           USE CC_Prop_Module
                           USE Iterative_Global
                           USE Pack_Global
                           IMPLICIT NONE
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          INTERFACE Cholesky
                   MODULE PROCEDURE Cholesky_d,                               &
                                    Cholesky_z
                          END INTERFACE Cholesky
!
                          INTERFACE Triangular_Solve
                   MODULE PROCEDURE Triangular_Solve_d,                       &
                                    Triangular_Solve_z
                          END INTERFACE Triangular_Solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           Contains
!***********************************************************************
!***********************************************************************
!deck Cholesky_d
!***begin prologue     Cholesky_d
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
!***end prologue       Cholesky_d
!
  SUBROUTINE Cholesky_d (input_matrix,output_matrix)
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)                :: input_matrix
  REAL*8, DIMENSION(:,:)                :: output_matrix
  INTEGER                               :: ifeig=.false.
!
  IF(print_parameter) THEN
     title='Input Matrix'
     call matprt(title,input_matrix,n3d,n3d,n3d,n3d,0,0,         &
                 rowlab,collab,1,eig,ifeig)
  END IF
  output_matrix(:,:) = input_matrix(:,:)
  call dpotrf('l',n3d,output_matrix,n3d,info)
  IF(print_parameter) THEN
     title='Cholesky Decomposition'
     call matprt(title,output_matrix,n3d,n3d,n3d,n3d,0,0,        &
                 rowlab,collab,1,eig,ifeig)
  END IF
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Cholesky_d
!***********************************************************************
!***********************************************************************
!deck Cholesky_z
!***begin prologue     Cholesky_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Cholesky decomposition of hermitian, positive
!***                   definite matrix.
!***references
!***routines called
!***end prologue       Cholesky_z
!
  SUBROUTINE Cholesky_z (input_matrix,output_matrix)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:,:)                :: input_matrix
  COMPLEX*16, DIMENSION(:,:)                :: output_matrix
!
  IF(print_parameter) THEN
     title='Input Matrix'
     call prntcmn(title,input_matrix,n3d,n3d,n3d,n3d,iout,'e')
  END IF
  output_matrix(:,:) = input_matrix(:,:)
  call zpotrf('l',n3d,output_matrix,n3d,info)
  IF(print_parameter) THEN
     title='Cholesky Decomposition'
     call prntcmn(title,output_matrix,n3d,n3d,n3d,n3d,iout,'e')
  END IF
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Cholesky_z
!***********************************************************************
!***********************************************************************
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
  SUBROUTINE Triangular_Solve_d (lower,solution,right_hand_side,packed)
  IMPLICIT NONE
  REAL*8,  DIMENSION(:,:)               :: lower
  REAL*8,  DIMENSION(:)                 :: solution
  REAL*8,  DIMENSION(:)                 :: right_hand_side
  INTEGER                               :: wptoin
  INTEGER                               :: first_row
  INTEGER                               :: last_row
  INTEGER                               :: ipnt
  LOGICAL                               :: packed
  INTEGER                               :: i, j, k 
  INTEGER                               :: count 
  INTEGER                               :: num 
  INTEGER                               :: ij, jk 
!
  ipnt = unit_pointer
  solution(:) = right_hand_side(:) 
  count = 0
  IF (packed) THEN
      IF(in_core) THEN
         DO i=1,n3d
            DO j=1,non_zero(i,ipnt)
               count = count + 1
               ij=row_buf(count,ipnt)
               solution(i) = solution(i)                                   &
                                         -                                 &
                             matrix_buf_d(count,ipnt) * solution(ij)
            END DO
            solution(i) = solution(i)/matrix_diag_d(i,ipnt)
         END DO
      ELSE
         Call iosys('rewind all on '//UNIT_NAME(ipnt)//'                   &
                     read-and-write',0,0,0,' ')
         Call iosys('read integer total_number_of_non_zero_elements from ' &
                     //UNIT_NAME(ipnt),1,num,0,' ')
         Call iosys('read integer non_zero_row_indices from '              &
                     //UNIT_NAME(ipnt),n3d,non_zero(1,ipnt),0,' ')
         Call iosys('read real diagonal_matrix_elements from '             &
                     //UNIT_NAME(ipnt),n3d,matrix_diag_d(1,ipnt),0,' ')
         Call iosys('read integer number_of_trips from '                   &
                     //UNIT_NAME(ipnt),1,trips,0,' ')
         Call iosys('read integer left_over from '                         &
                     //UNIT_NAME(ipnt),1,left_over,0,' ')
         first_row = 0
         last_row = 0
         DO i=1,trips
            first_row = last_row + 1
            last_row = last_row + max_row(ipnt)
            IF(i == trips) THEN
               last_row = last_row + left_over
            END IF
            Call iosys('read integer matrix_buffers from '                 &
                       //UNIT_NAME(ipnt)//' without rewinding',1,count,0,' ')
            Call iosys('read integer matrix_buffers from '                 &
                       //UNIT_NAME(ipnt)//' without rewinding',            &
                        count,row_buf(1,ipnt),0,' ')
            Call iosys('read integer matrix_buffers from '                 &
                       //UNIT_NAME(ipnt)//' without rewinding',            &
                        wptoin(count),matrix_buf_d(1,ipnt),0,' ')
            count = 0
            DO j=first_row, last_row
               DO k=1,non_zero(j,ipnt)
                  count = count + 1
                  jk=row_buf(count,ipnt)
                  solution(j) = solution(j)                                &
                                             -                             &
                                matrix_buf_d(count,ipnt) * solution(jk)
               END DO
               solution(j) = solution(j)/matrix_diag_d(j,ipnt)
            END DO
         END DO
      END IF
  ELSE
      DO i=1,n3d
         DO j=1,i-1
            solution(i) = solution(i)                                      &
                                       -                                   &
                          lower(i,j) * solution(j)
         END DO
         solution(i) = solution(i)/lower(i,i)
      END DO
  END IF
  IF(print_parameter) THEN
     title='Solution Vector'
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
  SUBROUTINE Triangular_Solve_z (lower,solution,right_hand_side,packed)
  IMPLICIT NONE
  COMPLEX*16,  DIMENSION(:,:)           :: lower
  COMPLEX*16,  DIMENSION(:)             :: solution
  COMPLEX*16,  DIMENSION(:)             :: right_hand_side
  INTEGER                               :: wptoin
  INTEGER                               :: first_row
  INTEGER                               :: last_row
  LOGICAL                               :: packed
  INTEGER                               :: ipnt
  INTEGER                               :: i, j, k 
  INTEGER                               :: count 
  INTEGER                               :: num 
  INTEGER                               :: ij, jk 
!
  ipnt = unit_pointer
  solution(:) = right_hand_side(:) 
  count = 0
  IF (packed) THEN
      IF(in_core) THEN
         DO i=1,n3d
            DO j=1,non_zero(i,ipnt)
               count = count + 1
               ij=row_buf(count,ipnt)
               solution(i) = solution(i)                                   &
                                         -                                 &
                             matrix_buf_z(count,ipnt) * solution(ij)
            END DO
            solution(i) = solution(i)/matrix_diag_z(i,ipnt)
         END DO
      ELSE
         Call iosys('rewind all on '//UNIT_NAME(ipnt)//'                   &
                     read-and-write',0,0,0,' ')
         Call iosys('read integer total_number_of_non_zero_elements from ' &
                     //UNIT_NAME(ipnt),1,num,0,' ')
         Call iosys('read integer non_zero_row_indices from '              &
                     //UNIT_NAME(ipnt),n3d,non_zero(1,ipnt),0,' ')
         Call iosys('read real diagonal_matrix_elements from '             &
                     //UNIT_NAME(ipnt),2*n3d,matrix_diag_z(1,ipnt),0,' ')
         Call iosys('read integer number_of_trips from '                   &
                     //UNIT_NAME(ipnt),1,trips,0,' ')
         Call iosys('read integer left_over from '                         &
                     //UNIT_NAME(ipnt),1,left_over,0,' ')
         first_row = 0
         last_row = 0
         DO i=1,trips
            first_row = last_row + 1
            last_row = last_row + max_row(ipnt)
            IF(i == trips) THEN
               last_row = last_row + left_over
            END IF
            Call iosys('read integer matrix_buffers from '                 &
                       //UNIT_NAME(ipnt)//' without rewinding',1,count,0,' ')
            Call iosys('read integer matrix_buffers from '                 &
                       //UNIT_NAME(ipnt)//' without rewinding',            &
                        count,row_buf(1,ipnt),0,' ')
            Call iosys('read integer matrix_buffers from '                 &
                       //UNIT_NAME(ipnt)//' without rewinding',            &
                        wptoin(2*count),matrix_buf_z(1,ipnt),0,' ')
            count = 0
            DO j=first_row, last_row
               DO k=1,non_zero(j,ipnt)
                  count = count + 1
                  jk=row_buf(count,ipnt)
                  solution(j) = solution(j)                                &
                                             -                             &
                                matrix_buf_z(count,ipnt) * solution(jk)
               END DO
               solution(j) = solution(j)/matrix_diag_z(j,ipnt)
            END DO
         END DO
      END IF
  ELSE
      DO i=1,n3d
         DO j=1,i-1
            solution(i) = solution(i)                                      &
                                       -                                   &
                          lower(i,j) * solution(j)
         END DO
         solution(i) = solution(i)/lower(i,i)
      END DO
  END IF
  IF(print_parameter) THEN
     title='Solution Vector'
     call prntfmn(title,solution,n3d,1,n3d,1,iout,'e')
  END IF
!***********************************************************************
!***********************************************************************
  END SUBROUTINE Triangular_Solve_z
!***********************************************************************
!***********************************************************************
  END  MODULE Preconditioner_Module
!***********************************************************************
!***********************************************************************
