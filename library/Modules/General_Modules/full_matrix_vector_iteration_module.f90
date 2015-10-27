!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! full_matrix_vector_iteration_module
!**begin prologue     full_matrix_vector_iteration_module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains the routines needed to multiply the
!***                  Hamiltonian on a vector or to mutiply the exponential
!***                  sector propagators on a vector. Explicit interfaces are 
!***                  used to allow a transparent use of generic subroutines 
!***                  which work for both real and complex vectors.  
!***                  This feature permits a single code to be used for both 
!***                  real and imaginary time propagation.
!***description       See the specific routined.
!**references
!**modules needed     See USE statements below
!**end prologue       full_matrix_vector_iteration_module
!***********************************************************************
!***********************************************************************
                      MODULE full_matrix_vector_iteration_module
                      USE input_output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      IMPLICIT NONE
  REAL*8                                        :: tol
  REAL*8                                        :: test_convergence
  INTEGER                                       :: n_iter
  INTEGER                                       :: the_size
  INTEGER                                       :: i
  INTEGER                                       :: j
  INTEGER                                       :: ii
  INTEGER                                       :: iter
  INTEGER                                       :: count
  LOGICAL                                       :: iterative_print
  CHARACTER(LEN=80)                             :: title
                      INTERFACE jacobi_iteration
              MODULE PROCEDURE jacobi_iteration_d,                     &
                               jacobi_iteration_z
                END INTERFACE  jacobi_iteration
!
                      INTERFACE gauss_seidel_iteration
              MODULE PROCEDURE gauss_seidel_iteration_d,               &
                               gauss_seidel_iteration_z
                END INTERFACE  gauss_seidel_iteration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck jacobi_iteration_d
!***begin prologue     jacobi_iteration_d    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           
!***author             schneider, b, i. (nsf)
!***source             
!***purpose            
!***                   
!
!***description        
!***                   
!------------------------------------------------------------------------------------
!
!                 !
!------------------------------------------------------------------------------------
!***references
!***routines called    
!
!***end prologue       jacobi_iteration_d                     
!
  SUBROUTINE jacobi_iteration_d (matrix,vector_in,vector_out,rhs)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                    :: matrix
  REAL*8, DIMENSION(:)                    :: vector_in
  REAL*8, DIMENSION(:)                    :: vector_out
  REAL*8, DIMENSION(:)                    :: rhs
  REAL*8                                  :: diff
  REAL*8                                  :: residual
!
!
!     Initial Guess
!
  count = 0
  DO i = 1, the_size
     count = count + 1
     vector_in(i) = rhs(i) / matrix(count)
     count = count + i
  END DO
!
!     Begin Iterations
!
  DO iter = 1,n_iter
     vector_out(:) = rhs(:)
!
!         Form V_o = b - M_od * V_i
!
     count = 0
     DO i =  1, the_size
        DO j = 1, i - 1
           count = count + 1
           vector_out(i) = vector_out(i) - matrix(count) * vector_in(j)
           vector_out(j) = vector_out(j) - matrix(count) * vector_in(i)
       END DO
       count = count + 1
     END DO
!         
!         Divide by Diagonal element to get new vector
!
     count = 0
     DO i = 1, the_size
        count = count + 1
        vector_out(i) = vector_out(i) / matrix(count)
        count = count + i
     END DO
     test_convergence = 0.d0
     DO i = 1, the_size
        diff = vector_out(i) - vector_in(i)
        test_convergence = test_convergence + diff * diff
     END DO
     test_convergence = sqrt(test_convergence)
     IF(iterative_print) THEN     
        Write(iout,1) iter, test_convergence
        title='Vector'
        Write(iout,2) title
        Write(iout,3) vector_out
     END IF
     IF (test_convergence <= tol ) THEN
         EXIT
     ELSE
        vector_in(:) = vector_out(:)
     END IF
  END DO
1 Format(/,5x,'Iteration = ',i5, 2x, 'Convergence = ',e15.8)
2 Format(a80)
3 Format(5e15.8)
!
 END SUBROUTINE jacobi_iteration_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!deck jacobi_iteration_z
!***begin prologue     jacobi_iteration_z    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           
!***author             schneider, b, i. (nsf)
!***source             
!***purpose            
!***                   
!
!***description        
!***                   
!------------------------------------------------------------------------------------
!
!                 !
!------------------------------------------------------------------------------------
!***references
!***routines called    
!
!***end prologue       jacobi_iteration_z                     
!
  SUBROUTINE jacobi_iteration_z (matrix,vector_in,vector_out,rhs)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)                :: matrix
  COMPLEX*16, DIMENSION(:)                :: vector_in
  COMPLEX*16, DIMENSION(:)                :: vector_out
  COMPLEX*16, DIMENSION(:)                :: rhs
  COMPLEX*16                              :: diff
!
  count = 0
  DO i = 1, the_size
     count = count + 1
     vector_in(i) = rhs(i) / matrix(count)
     count = count + i
  END DO
  DO iter = 1,n_iter
     vector_out(:) = rhs(:)
     count = 0
     DO i =  1, the_size
        DO j = 1, i - 1
           count = count + 1
           vector_out(i) = vector_out(i) - matrix(count) * vector_in(j)
           vector_out(j) = vector_out(j) - conjg(matrix(count)) * vector_in(i)
       END DO
       count = count + 1
     END DO
     count = 0
     DO i = 1, the_size
        count = count + 1
        vector_out(i) = vector_out(i) / matrix(count)
        count = count + i
     END DO
     test_convergence = 0.d0
     DO i = 1, the_size
        diff = vector_out(i) - vector_in(i)
        test_convergence = test_convergence + diff * conjg(diff)
     END DO
     test_convergence = sqrt(test_convergence)
     IF(iterative_print) THEN     
        Write(iout,1) iter, test_convergence
        title='Vector'
        Write(iout,2) title
        Write(iout,3) vector_out
     END IF
     IF (test_convergence <= tol ) THEN
         EXIT
     ELSE
        vector_in(:) = vector_out(:)
     END IF
  END DO
1 Format(/,5x,'Iteration = ',i5, 2x, 'Convergence = ',e15.8)
2 Format(a80)
3 Format(5e15.8)
!
 END SUBROUTINE jacobi_iteration_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck gauss_seidel_iteration_d
!***begin prologue     gauss_seidel_iteration_d    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           
!***author             schneider, b, i. (nsf)
!***source             
!***purpose            
!***                   
!
!***description        
!***                   
!------------------------------------------------------------------------------------
!
!                 !
!------------------------------------------------------------------------------------
!***references
!***routines called    
!
!***end prologue       gauss_seidel_iteration_d                     
!
  SUBROUTINE gauss_seidel_iteration_d (matrix,vector_in,vector_out,rhs,omega)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                    :: matrix
  REAL*8, DIMENSION(:)                    :: vector_in
  REAL*8, DIMENSION(:)                    :: vector_out
  REAL*8, DIMENSION(:)                    :: rhs
  REAL*8                                  :: diff
  REAL*8                                  :: omega
  REAL*8                                  :: x_i
  count = 0
  DO i = 1, the_size
     count = count + 1
     vector_out(i) = rhs(i) / matrix(count)
     count = count + i
  END DO
!
  vector_in(1:the_size) = vector_out(1:the_size)
  DO iter = 1,n_iter
     DO i =  1, the_size
        x_i = rhs(i)
        ii = i*(i-1)/2 
        DO j = 1, i-1
           count = ii +j 
           x_i = x_i - matrix(count) * vector_out(j)
        END DO
        DO j = i+1, the_size
           count = j*(j-1)/2 + i
           x_i = x_i - matrix(count) * vector_out(j)
        END DO
        ii = ii + i
        vector_out(i) = (1.d0 - omega) * vector_out(i) + omega * x_i / matrix(ii)
     END DO
     test_convergence = 0.d0
     DO i = 1, the_size
        diff = vector_out(i) - vector_in(i)
        test_convergence = test_convergence + diff * diff
     END DO
     test_convergence = sqrt(test_convergence)
     IF(iterative_print) THEN     
        Write(iout,1) iter, test_convergence
        title='Vector'
        Write(iout,2) title
        Write(iout,3) vector_out
     END IF
     IF (test_convergence <= tol ) THEN
         Write(iout,4) iter, test_convergence
         EXIT
     ELSE
        vector_in(1:the_size) = vector_out(1:the_size)
     END IF
  END DO
  vector_in(:) = 0.d0
  count = 0
  DO i=1,the_size
     DO j=1,i-1
        count = count + 1
        vector_in(i) = vector_in(i) + matrix(count) * vector_out(j)
        vector_in(j) = vector_in(j) + matrix(count) * vector_out(i)
     END DO
     count = count + 1
     vector_in(i) = vector_in(i) + matrix(count) * vector_out(i)     
  END DO
  vector_in(1:the_size) = vector_in(1:the_size) - rhs(1:the_size)
  write(iout,*) vector_in
1 Format(/,5x,'Iteration = ',i5, 2x, 'Convergence = ',e15.8)
2 Format(a80)
3 Format(5e15.8)
4 Format(/,5x,'Iteration = ',i4,2x,'Convergence = ',e15.8)
!
!
 END SUBROUTINE gauss_seidel_iteration_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!deck gauss_seidel_iteration_z
!***begin prologue     gauss_seidel_iteration_z    
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           
!***author             schneider, b, i. (nsf)
!***source             
!***purpose            
!***                   
!
!***description        
!***                   
!------------------------------------------------------------------------------------
!
!                 !
!------------------------------------------------------------------------------------
!***references
!***routines called    
!
!***end prologue       gauss_seidel_iteration_z                     
!
  SUBROUTINE gauss_seidel_iteration_z (matrix,vector_in,vector_out,rhs,omega)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)                :: matrix
  COMPLEX*16, DIMENSION(:)                :: vector_in
  COMPLEX*16, DIMENSION(:)                :: vector_out
  COMPLEX*16, DIMENSION(:)                :: rhs
  COMPLEX*16                              :: diff
  REAL*8                                  :: omega
!
  count = 0
  DO i = 1, the_size
     count = count + 1
     vector_in(i) = rhs(i) / matrix(count)
     count = count + i
  END DO
  DO iter = 1,n_iter
     vector_out(:) = vector_in(:)
     DO i =  1, the_size
        ii = i*(i-1)/2 
        DO j = 1, i-1
           count = ii +j 
           vector_out(i) = rhs(i) - conjg(matrix(count)) * vector_out(j)
        END DO
        DO j = i+1, the_size
           count =  i + j*(j-1)/2
           vector_out(i) = rhs(i) - matrix(count) * vector_in(j)
        END DO
     END DO
     count = 0
     DO i = 1, the_size
        count = count + 1
        vector_out(i) = vector_out(i) / matrix(count)
        count = count + i
     END DO
     test_convergence = 0.d0
     DO i = 1, the_size
        diff = vector_out(i) - vector_in(i)
        test_convergence = test_convergence + diff * conjg(diff)
     END DO
     test_convergence = sqrt(test_convergence)
     IF(iterative_print) THEN     
        Write(iout,1) iter, test_convergence
        title='Vector'
        Write(iout,2) title
        Write(iout,3) vector_out
     END IF
     IF (test_convergence <= tol ) THEN
         EXIT
     ELSE
        vector_in(:) = vector_out(:)
     END IF
  END DO
!
!
1 Format(/,5x,'Iteration = ',i5, 2x, 'Convergence = ',e15.8)
2 Format(a80)
3 Format(5e15.8)
 END SUBROUTINE gauss_seidel_iteration_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END       MODULE full_matrix_vector_iteration_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
