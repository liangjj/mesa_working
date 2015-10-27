!***********************************************************************
! Read_Matrix_From_Input
!**begin prologue     Read_Matrix_From_Input
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Davidson, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains all of the major subroutines to 
!***                  solve a large eigenvalue problem or a large set of
!***                  linear equations using the Davidson algorithm.  
!***                  Explicit interfaces are used to allow
!***                  a transparent use of generic subroutines which work
!***                  for both real and complex vectors.  This feature
!***                  permits a single code to be used for both real symmetric
!***                  and Hermitian matrices.
!***description       Given a starting vector, a number of iterations
!***                  are performed until the desired eigenvalues or linear
!***                  system has been solved to a given accuracy criterion.  
!**references
!**modules needed     See USE statements below
!**comments           In this portable version I have disabled all unnecessary
!**                   writing to files.  The original Fortran is commented out.
!**                   In addition, there is no option to compute the autocorrelation
!**                   function as this would require reading and manipulating the
!**                   initial state wavefunction from a file.
!**end prologue       Read_Matrix_From_Input
!***********************************************************************
!***********************************************************************
  SUBROUTINE Read_Matrix_From_Input
                           USE Davidson_Module 
                           USE Iterative_Global 
  IMPLICIT NONE
  INTEGER                                        :: i
  INTEGER                                        :: j
  LOGICAL                                        :: dollar
  CHARACTER(LEN=3)                               :: itoc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  maximum_number_of_non_zero_elements = matrix_size * ( matrix_size - 1) / 2
  ALLOCATE(ibuf(2,maximum_number_of_non_zero_elements))
  IF(matrix_type == 'real_symmetric') THEN
     ALLOCATE(matrix_d(matrix_size,matrix_size), h_buf_d(maximum_number_of_non_zero_elements), diag_d(matrix_size))
     IF (type_calculation == 'linear_system') THEN
         ALLOCATE(rhs_d(matrix_size))
     END IF
  ELSE IF(matrix_type == 'hermitian') THEN
     ALLOCATE(matrix_z(matrix_size,matrix_size), h_buf_z(maximum_number_of_non_zero_elements), diag_d(matrix_size))
     IF (type_calculation == 'linear_system') THEN
         ALLOCATE(rhs_z(matrix_size))
     END IF
  END IF
  IF( dollar ('$matrix_input',card,cpass,inp) ) THEN
      IF(matrix_type == 'real_symmetric') THEN      
         DO i=1,matrix_size
            call fparr(card,'matrix_row_'//itoc(i),h_buf_d,i,' ')
            matrix_d(i,1:i) = h_buf_d(1:i)
         END DO
         DO i = 1, matrix_size
            DO j = 1, i
               matrix_d(j,i) = matrix_d(i,j)
            END DO
         END DO
         title='input matrix'
         Call prntfmn(title,matrix_d,matrix_size,matrix_size,matrix_size,matrix_size,iout,'e')
         IF (type_calculation == 'linear_system') THEN
             call fparr(card,'right_hand_side',rhs_d,matrix_size,' ')
             title='input rhs'
             Call prntfmn(title,rhs_d,matrix_size,1,matrix_size,1,iout,'e')
         END IF
         Call Pack_Matrix(matrix_d)
         DEALLOCATE(matrix_d)
      ELSE IF(matrix_type == 'hermititan') THEN      
         DO i=1,matrix_size
            call fparr(card,'matrix_row_'//itoc(i),h_buf_z,i+i,' ')
            matrix_z(i,1:i) = h_buf_z(1:i)
         END DO
         DO i = 1, matrix_size
            DO j = 1, i
               matrix_z(j,i) = conjg (matrix_z(i,j) )
            END DO
         END DO
         title='input matrix'
         Call prntcmn(title,matrix_z,matrix_size,matrix_size,matrix_size,matrix_size,iout,'e')
         IF (type_calculation == 'linear_system') THEN
             call fparr(card,'right_hand_side',rhs_z,matrix_size+matrix_size,' ')
             title='input rhs'
             Call prntcmn(title,rhs_z,matrix_size,1,matrix_size,1,iout,'e')
         END IF
         Call Pack_Matrix(matrix_z)
         DEALLOCATE(matrix_z)
      END IF
  END IF
  END SUBROUTINE Read_Matrix_From_Input
!***********************************************************************
!***********************************************************************




