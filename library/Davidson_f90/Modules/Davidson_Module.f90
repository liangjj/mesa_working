!***********************************************************************
! Davidson_Module
!**begin prologue     Davidson_Module
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
!**end prologue       Davidson_Module
!***********************************************************************
!***********************************************************************
                           MODULE Davidson_Module
                           USE Iterative_Global
                           USE packed_matrix_vector_multiply_module
                      IMPLICIT NONE
  INTEGER                                :: matrix_size
  INTEGER                                :: number_of_roots
  INTEGER                                :: number_of_right_hand_sides
  INTEGER                                :: number_of_passes
  INTEGER                                :: number_of_roots_per_pass
  INTEGER                                :: maximum_number_of_iterations
  INTEGER                                :: maximum_number_of_non_zero_elements
  INTEGER                                :: maximum_number_of_davidson_vectors
  INTEGER                                :: guess_size
  INTEGER                                :: number_of_guess_vectors
  INTEGER                                :: n_tri
  INTEGER                                :: number_converged
  INTEGER                                :: number_unconverged
  REAL, DIMENSION(20)                    :: c_time
  REAL*8                                 :: davidson_convergence
  REAL*8                                 :: overlap_tolerance
  REAL*8                                 :: dvd_error
  REAL*8                                 :: maximum_error
  LOGICAL, DIMENSION(20)                 :: prdvd
  LOGICAL                                :: dvdall
  INTEGER                                :: total_number_of_iterations
  CHARACTER (LEN=16)                     :: matrix_type
  CHARACTER (LEN=16)                     :: matrix_source
  CHARACTER (LEN=80)                     :: preconditioner
  CHARACTER (LEN=80)                     :: storage_mode
  CHARACTER (LEN=8)                      :: convergence_control
  REAL*8, DIMENSION(:), ALLOCATABLE      :: b_d
  COMPLEX*16, DIMENSION(:), ALLOCATABLE  :: b_z

!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                           INTERFACE Allocation
             MODULE PROCEDURE ALlocation
                       END INTERFACE Allocation
!
                           INTERFACE Pack_matrix
             MODULE PROCEDURE Pack_matrix_d,                       &
                              Pack_matrix_z
                       END INTERFACE Pack_matrix
!
                           INTERFACE Gram_Schmidt
             MODULE PROCEDURE Gram_Schmidt_d,                      &
                              Gram_Schmidt_z
                       END INTERFACE Gram_Schmidt
!
                           INTERFACE Initialize
             MODULE PROCEDURE Initialize_d,                        &
                              Initialize_z
                       END INTERFACE Initialize
!
                           INTERFACE Linear_System_Driver
             MODULE PROCEDURE Linear_System_Driver_d,              &
                              Linear_System_Driver_z
                       END INTERFACE Linear_System_Driver
!
                           INTERFACE Solution_Vector
             MODULE PROCEDURE Solution_Vector_d,                   &
                              Solution_Vector_z
                       END INTERFACE Solution_Vector
!
                           INTERFACE Solve_Small_Linear_System
             MODULE PROCEDURE Solve_Small_Linear_System_d,          &              
                              Solve_Small_Linear_System_z              
                       END INTERFACE Solve_Small_Linear_System
!
                           INTERFACE Solve_Small_Eigenvalue_Problem
             MODULE PROCEDURE Solve_Small_Eigenvalue_Problem_d,     &              
                              Solve_Small_Eigenvalue_Problem_z                            
                       END INTERFACE Solve_Small_Eigenvalue_Problem                             
!
                           INTERFACE Convergence_Test
             MODULE PROCEDURE Convergence_Test_d,                  &
                              Convergence_Test_z  
                       END INTERFACE Convergence_Test
!
                           INTERFACE Guess_Solution
             MODULE PROCEDURE Guess_Solution_d,                    &
                              Guess_Solution_z  
                       END INTERFACE Guess_Solution
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Davidson_Data
!***begin prologue     Davidson_Data
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            data entry for davidson routine.
!***
!***references

!***routines called
!***end prologue       dvdddat

  SUBROUTINE Davidson_Data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LOGICAL                   :: dollar
  LOGICAL                   :: logkey
  CHARACTER(LEN=80)         :: chrkey
  INTEGER                   :: intkey
  REAL*8                    :: fpkey
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF( dollar ('$davidson_parameters',card,cpass,inp) ) THEN
      overlap_tolerance = fpkey(card,'overlap_tolerance',1.0D-10,' ')
      matrix_size = intkey(card,'matrix_size',2,' ')
      number_of_roots = intkey(card,'number_of_roots',1,' ')
      number_of_roots = intkey(card,'number_of_right_hand_sides',1,' ')
      number_of_roots_per_pass = intkey(card,'number_of_roots_per_pass',     &
                                 number_of_roots,' ')
      maximum_number_of_iterations = intkey(card,'maximum_number_of_iterations',1,' ')
      maximum_number_of_davidson_vectors = intkey(card,                      &
                                           'maximum_number_of_davidson_vectors',1,' ')
      guess_size = intkey(card,'guess_size',10*number_of_roots_per_pass,' ')
      guess_size = min(matrix_size,guess_size)
      davidson_convergence = fpkey(card,'davidson_convergence',1.0D-08,' ')
      drop_tol = fpkey(card,'drop_tolerance',1.0d-10,' ')
      number_of_guess_vectors = intkey(card,'number_of_guess_vectors',guess_size,' ')
      number_of_guess_vectors = min(1,guess_size,maximum_number_of_davidson_vectors/2)
      matrix_type = chrkey(card,'matrix_type','real_symmetric',' ')
      matrix_source = chrkey(card,'matrix_source','from_input',' ')
      type_calculation = chrkey(card,'type_calculation','linear_system',' ')
      preconditioner = chrkey(card,'preconditioner','diagonal',' ')
  END IF
  prdvd(1)=logkey(card,'print=trials',.false.,' ')
  prdvd(2)=logkey(card,'print=vectors',.false.,' ')
  prdvd(3)=logkey(card,'print=h_on_vectors',.false.,' ')
  prdvd(4)=logkey(card,'print=hamiltonian',.false.,' ')
  prdvd(5)=logkey(card,'print=iteration_information',.false.,' ')
  prdvd(6)=logkey(card,'print=residuals',.false.,' ')
  prdvd(7)=logkey(card,'print=transformed_vectors',.false.,' ')
  prdvd(8)=logkey(card,'print=transformed_h_on_vectors', .false.,' ')
  prdvd(9)=logkey(card,'print=new_trial_vectors',.false.,' ')
  prdvd(10)=logkey(card,'print=overlaps',.false.,' ')
  prdvd(11)=logkey(card,'print=guess_matrix',.false.,' ')
  prdvd(12)=logkey(card,'print=guess_vectors',.false.,' ')
  prdvd(13)=logkey(card,'print=packed_matrix',.false.,' ')
  prdvd(14)=logkey(card,'print=residuals',.false.,' ')
  dvdall=logkey(card,'print=all',.false.,' ')
  prdvd(:)=.false.
  IF(dvdall) THEN
     prdvd(:)=.true.
  END IF
  WRITE(iout,1) matrix_size, guess_size, number_of_guess_vectors, maximum_number_of_davidson_vectors,        &
                maximum_number_of_iterations, overlap_tolerance, davidson_convergence
  IF (type_calculation == 'eigenvalue') THEN
      WRITE(iout,2) number_of_roots, number_of_roots_per_pass
  END IF
1 FORMAT(/,15X,'iterative diagonalization information',/,/,5X,      &
               'size of matrix                      = ',i8,/,5X,    &
               'size of guess matrix                = ',i7,/,5X,    &
               'size of initial vectors             = ',i7,/,5X,    &
               'maximum number of davidson vectors  = ',i7,/,5X,    &
               'max. number of iterations           = ',i7,/,5X,    &
               'overlap tolerance                   = ',e15.8,/,5X, &
               'convergence criterion               = ',e15.8)
2 FORMAT(/,15X,'number of roots                     = ',i4,/,5X,  &
               'number of roots at a time           = ',i4)


  END SUBROUTINE Davidson_Data
!***********************************************************************
!***********************************************************************
!deck Allocation
!***begin prologue     Allocation
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            allocation
!***
!***references

!***routines called
!***end prologue       Allocation

  SUBROUTINE Allocation
  INTEGER                                 :: i_plus  
!
! At this point we have allocated and computed the arrays needed for the
! packed Hamiltonian, the guess eigenvectors or guess to the linear system
! and the right hand side.
  i_plus = maximum_number_of_davidson_vectors + 1
  n_tri = i_plus * ( i_plus + 1 ) / 2
  Write(iout,1)
  IF (matrix_type =='real_symmetric') THEN
      IF ( type_calculation == 'linear_system') THEN
          ALLOCATE(v_in_d(matrix_size),                                           &
                   v_out_d(matrix_size),                                          &
                   vec_d(matrix_size,0:maximum_number_of_davidson_vectors),       &
                   h_vectors_d(matrix_size,0:maximum_number_of_davidson_vectors), &
                   h_mat_tri_d(0:n_tri), h_mat_work_tri_d(0:n_tri),               &
                   small_rhs_d(0:maximum_number_of_davidson_vectors),             & 
                   small_rhs_work_d(0:maximum_number_of_davidson_vectors),        & 
                   ipvt(i_plus))
      ELSE IF (type_calculation == 'eigenvalues') THEN
          ALLOCATE(vec_d(matrix_size,0:maximum_number_of_davidson_vectors),       &
                   h_vectors_d(matrix_size,0:maximum_number_of_davidson_vectors), &
                   h_mat_d(0:maximum_number_of_davidson_vectors,                  &
                           0:maximum_number_of_davidson_vectors),                 &
                   h_mat_work_d(0:maximum_number_of_davidson_vectors,             &
                                0:maximum_number_of_davidson_vectors),            &
                   eigen_values(0:maximum_number_of_davidson_vectors),            & 
                   working_eigen_values(0:maximum_number_of_davidson_vectors),    & 
                   work_d(5*i_plus) )
      END IF
  ELSE IF(matrix_type =='hermitian') THEN
      IF ( type_calculation == 'linear_system') THEN
          ALLOCATE(v_in_z(matrix_size),                                           &
                   v_out_z(matrix_size),                                          &
                   vec_z(matrix_size,0:maximum_number_of_davidson_vectors),       &
                   h_vectors_z(matrix_size,0:maximum_number_of_davidson_vectors), &
                   h_mat_tri_z(0:n_tri), h_mat_work_tri_z(0:n_tri),               &
                   small_rhs_z(0:maximum_number_of_davidson_vectors),             & 
                   small_rhs_work_z(0:maximum_number_of_davidson_vectors),        & 
                   ipvt(i_plus))
      ELSE IF (type_calculation == 'eigenvalues') THEN
          ALLOCATE(vec_z(matrix_size,0:maximum_number_of_davidson_vectors),       &
                   h_vectors_z(matrix_size,0:maximum_number_of_davidson_vectors), &
                   h_mat_z(0:maximum_number_of_davidson_vectors,                  &
                           0:maximum_number_of_davidson_vectors),                 &
                   h_mat_work_z(0:maximum_number_of_davidson_vectors,             &
                                0:maximum_number_of_davidson_vectors),            &
                   eigen_values(0:maximum_number_of_davidson_vectors),            & 
                   working_eigen_values(0:maximum_number_of_davidson_vectors),    & 
                   work_z(5*i_plus), rwork(5*i_plus) )
      END IF
  END IF
1 Format(/,10x,'Allocating Main Computational Arrays')
  END SUBROUTINE Allocation
!***********************************************************************
!***********************************************************************
!deck Pack_Matrix_d
!***begin prologue     Pack_Matrix_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Pack the Matrix.
!***
!***references

!***routines called
!***end prologue       Pack_Matrix_d
  SUBROUTINE Pack_Matrix_d(matrix)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL*8, DIMENSION(:,:)                :: matrix
  INTEGER                               :: i 
  INTEGER                               :: j 
!
!     pack all non-zero, non-diagonal elements in the lower triangle
!
  Write(iout,1)
  non_zero=0
  DO i=1,matrix_size
     DO j=1,i-1
        IF(abs(matrix(i,j)) > drop_tol) THEN
           non_zero=non_zero + 1
           ibuf(1,non_zero) = i
           ibuf(2,non_zero) = j
           h_buf_d(non_zero)=matrix(i,j)
        END IF
     END DO
  END DO 
!
!
  DO i=1,matrix_size
     diag_d(i)= matrix(i,i)
  END DO
  Write(iout,2) non_zero
  IF (prdvd(13)) THEN
      DO i=1,non_zero
         write(iout,3) ibuf(1,i), ibuf(2,i), h_buf_d(i)
      END DO
      Write(iout,4) 
      Write(iout,5) diag_d(1:matrix_size)
  END IF
1 Format(/,5x,'Packing the Matrix')  
2 Format(/,5x,'Number Non Zero Elements = ', i5)  
3 Format(5x,'I = ',i5,2x,'J = ',i5,2x,'Element = ',e15.8)
4 Format(/,5x,'Diagonal')
5 Format(5e15.8)
  END SUBROUTINE Pack_Matrix_d
!***********************************************************************
!***********************************************************************
!deck Pack_Matrix_z
!***begin prologue     Pack_Matrix_z
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Pack the Matrix.
!***
!***references

!***routines called
!***end prologue       Pack_Matrix_z
  SUBROUTINE Pack_Matrix_z(matrix)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  COMPLEX*16, DIMENSION(:,:)            :: matrix
  INTEGER                               :: i 
  INTEGER                               :: j 
!
!     pack the complex conjugate of all non-zero, non-diagonal elements in the lower triangle
!
  Write(iout,1)
  non_zero=0
  DO i=1,matrix_size
     DO j=1,i-1
        IF(abs(matrix(i,j)).gt.drop_tol) THEN
           non_zero=non_zero + 1
           ibuf(1,non_zero) = i
           ibuf(2,non_zero) = j
           h_buf_z(non_zero) = matrix(i,j)
        END IF
     END DO
  END DO 
!
! The matrix is assumed to be Hermitian.  Thus the diagonal elements must be real.  It is possible that they
! have some small imaginary part which is here removed using the real part of the variable.
!
  DO i=1,matrix_size
     diag_d(i) = Real( matrix(i,i) )
  END DO
  Write(iout,2) non_zero
  IF (prdvd(13)) THEN
      DO i=1,non_zero
         write(iout,3) ibuf(1,i), ibuf(2,i), h_buf_z(i)
      END DO
      Write(iout,4) 
      Write(iout,5) diag_d(1:matrix_size)
  END IF
1 Format(/,5x,'Packing the Matrix')  
2 Format(/,5x,'Number Non Zero Elements = ', i5)  
3 Format(5x,'I = ',i5,2x,'J = ',i5,2x,'Element = ',2e15.8)
4 Format(/,5x,'Diagonal')
5 Format(5e15.8)
  END SUBROUTINE Pack_Matrix_z
!***********************************************************************
!***********************************************************************
!deck @(#)Guess_Solution_d
  SUBROUTINE Guess_Solution_d(guess_vectors,guess_rhs)
!***begin prologue     Guess_Solution_d
!***date written       0803083   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           guess
!***author             
!***source             Guess_Solution
!***purpose            to form a portion of the matrix and extract starting solutions
!***description
!***references
!***routines called    (none)
!***end prologue       Guess_Solution_d
  IMPLICIT NONE
  REAL*8, DIMENSION(matrix_size,guess_size)   :: guess_vectors
  REAL*8, DIMENSION(matrix_size)              :: guess_rhs
  INTEGER                                     :: i
  INTEGER                                     :: j
  INTEGER                                     :: ii
  INTEGER                                     :: jj
  INTEGER                                     :: ref_walk
  INTEGER                                     :: i_guess
  INTEGER                                     :: j_guess
  INTEGER                                     :: i_tri
  INTEGER                                     :: count
  REAL*8, DIMENSION(:,:), ALLOCATABLE         :: guess_matrix
  REAL*8, DIMENSION(:), ALLOCATABLE           :: eigv
  REAL*8, DIMENSION(:), ALLOCATABLE           :: temp
  REAL*8, DIMENSION(:), ALLOCATABLE           :: work
  REAL*8                                      :: e_guess
  REAL*8                                      :: tmp
  INTEGER, DIMENSION(:), ALLOCATABLE          :: guess_pointers
  INTEGER, DIMENSION(:), ALLOCATABLE          :: ipivot
!
  write(iout,1) guess_size
! This is an attempt to pick out the part of the matrix with the largest diagonals
! It will not always be best, but it cannot be worse.
!
  ALLOCATE(temp(matrix_size), guess_pointers(matrix_size))
  temp(1:matrix_size) = diag_d(1:matrix_size)
  guess_pointers(:) = 0
  DO i = 1, guess_size
     e_guess = 1.d+30
     DO j = 1, matrix_size
        IF(temp(j) <  e_guess) THEN
           ref_walk = j
           e_guess = temp(j)
        END IF
     END DO
     temp(ref_walk) = 1.0d+30
     guess_pointers(ref_walk) = i
  END DO
!
  IF ( type_calculation == 'linear_system') THEN
       i_tri = guess_size * ( guess_size + 1 ) / 2       
       ALLOCATE( h_mat_tri_d(i_tri) )
       h_mat_tri_d = 0.d0
       DO i = 1, non_zero
          i_guess = guess_pointers(ibuf(1,i))
          IF (i_guess <= 0 ) cycle
          j_guess = guess_pointers(ibuf(2,i))
          IF (j_guess <= 0 ) cycle
          ii = max(i_guess,j_guess)
          jj = min(i_guess,j_guess)
          i_tri = ii * ( ii - 1 ) / 2 + jj      
          h_mat_tri_d(i_tri) = h_mat_tri_d(i_tri)  + h_buf_d(i)
       END DO
       DO  i=1,matrix_size
           i_guess = guess_pointers(i)
           IF (i_guess <= 0) cycle
           i_tri = i_guess * ( i_guess + 1 ) / 2
           h_mat_tri_d(i_tri) = diag_d(i)
       END DO
       IF( prdvd(11) ) THEN
           i_tri = 0
           title='Guess Matrix'
           write(iout,2) title
           DO i=1,guess_size
              write(iout,3) i
              write(iout,4) (h_mat_tri_d(j), j=i_tri + 1, i_tri + i )
              i_tri = i_tri + i
           END DO
       END IF
       ALLOCATE(ipivot(guess_size))
       DO  i=1,matrix_size 
           i_guess = guess_pointers(i)
           IF (i_guess <= 0) cycle
           guess_rhs(i_guess) = rhs_d(i)
       END DO
       call dspsv('u',guess_size,1,h_mat_tri_d,ipivot,guess_rhs,guess_size,info)
       temp(1:matrix_size) = 0.d0       
       tmp = 0.0D+00
       DO  i = 1, matrix_size
           IF (guess_pointers(i) > 0) THEN
               temp(i) = guess_rhs(guess_pointers(i))
               IF (ABS(temp(i)) > tmp) THEN
                   tmp = ABS(temp(i))
               END IF
           END IF
       END DO
       guess_rhs(1:matrix_size) = temp(1:matrix_size)
       Call iosys('open guess_rhs as new',0,0,0,'Guess_RHS.out')
       Call iosys('write real Guess_RHS to Guess_RHS',matrix_size,guess_rhs,0,' ')
       Call iosys('read real Guess_RHS from Guess_RHS',matrix_size,guess_rhs,0,' ')
       IF (prdvd(12)) THEN
           title = 'Guess Solution'
           CALL prntfmn(title,guess_rhs,matrix_size,1,matrix_size,1,iout,'e')
       END IF
       DEALLOCATE(temp, guess_pointers, h_mat_tri_d, ipivot)
  ELSE IF( type_calculation == 'eigenvalues') THEN
       Call iosys('open guess_vectors as new',0,0,0,'Trial_Vectors.out')
       ALLOCATE( h_mat_d(guess_size,guess_size) )
       h_mat_d = 0.d0
       DO i = 1, non_zero
          i_guess = guess_pointers(ibuf(1,i))
          IF (i_guess <= 0 ) cycle
          j_guess = guess_pointers(ibuf(2,i))
          IF (j_guess <= 0 ) cycle
          ii = max(i_guess,j_guess)
          jj = min(i_guess,j_guess)
          h_mat_d(ii,jj) = h_mat_d(ii,jj) + h_buf_d(i)
       END DO
       DO  i=1,matrix_size
           i_guess = guess_pointers(i)
           IF (i_guess <= 0) cycle
           h_mat_d(i_guess,i_guess) = diag_d(i)
       END DO
       DO i = 1, guess_size
          DO j = 1, i
             h_mat_d(j,i) = h_mat_d(i,j)
          END DO
       END DO
       IF (prdvd(11)) THEN
           title = 'Guess Matrix'
           CALL prntfmn(title,h_mat_d,guess_size,guess_size,guess_size,guess_size,iout,'e')
        END IF
       ALLOCATE( work(5*guess_size), eigv(guess_size))
       CALL dsyev('v','l',guess_size,h_mat_d,guess_size,eigv,work,5*guess_size,info)
       DO  i = 1, number_of_guess_vectors
           temp(1:matrix_size) = 0.d0       
           tmp = 0.0D+00
           DO  j = 1, matrix_size
               IF (guess_pointers(j) > 0) THEN
                   temp(j) = h_mat_d(guess_pointers(j),i)
                   IF (ABS(temp(j)) > tmp) THEN
                       tmp = ABS(temp(j))
                   END IF
               END IF
           END DO
           guess_vectors(1:matrix_size,i) = temp(1:matrix_size)
       END DO
       IF (prdvd(12)) THEN
           title = 'Guess Eigenvalues'
           CALL prntfmn(title,eigv,guess_size,1,guess_size,1,iout,'e')
           title = 'Guess Eigenvectors'
           CALL prntfmn(title,h_mat_d,guess_size,guess_size,guess_size,guess_size,iout,'e')
       END IF
       Call iosys('write real guess_eigenvalues to guess_vectors',guess_size,eigv,0,' ')
       Call iosys('write real guess_eigenvectors to guess_vectors',guess_size*matrix_size,h_mat_d,0,' ')
       DEALLOCATE(temp, guess_pointers, h_mat_d, eigv, work)
  END IF
1    FORMAT (/,t15,'Results from Guess Matrix = ', i5)
2 Format(a80)
3 Format(/,5x,'Row = ',i4)
4 Format( (15x,5f10.5) )
END SUBROUTINE Guess_Solution_d
!***********************************************************************
!***********************************************************************
!deck @(#)Guess_Solution_z
  SUBROUTINE Guess_Solution_z(guess_vectors,guess_rhs)
!***begin prologue     Guess_Solution_z
!***date written       0803083   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           guess
!***author             
!***source             Guess_Solution
!***purpose            to form a portion of the matrix and extract starting solutions
!***description
!***references
!***routines called    (none)
!***end prologue       Guess_Solution_z
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(matrix_size,guess_size)   :: guess_vectors
  COMPLEX*16, DIMENSION(matrix_size)              :: guess_rhs
  INTEGER                                         :: i
  INTEGER                                         :: j
  INTEGER                                         :: ii
  INTEGER                                         :: jj
  INTEGER                                         :: ref_walk
  INTEGER                                         :: i_guess
  INTEGER                                         :: j_guess
  INTEGER                                         :: i_tri
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE         :: guess_matrix
  REAL*8, DIMENSION(:), ALLOCATABLE               :: eigv
  REAL*8, DIMENSION(:), ALLOCATABLE               :: temp
  COMPLEX*16, DIMENSION(:), ALLOCATABLE           :: work
  REAL*8, DIMENSION(:), ALLOCATABLE               :: rwork
  REAL*8                                          :: e_guess
  REAL*8                                          :: tmp
  INTEGER, DIMENSION(:), ALLOCATABLE              :: guess_pointers
  INTEGER, DIMENSION(:), ALLOCATABLE              :: ipivot
!
  write(iout,1) guess_size
! This is an attempt to pick out the part of the matrix with the largest diagonals
! It will not always be best, but it cannot be worse.
!
  ALLOCATE(temp(matrix_size), guess_pointers(matrix_size))
  temp(1:matrix_size) = diag_d(1:matrix_size) 
  DO i = 1, guess_size
     e_guess = 1.d+30
     DO j = 1, matrix_size
        IF(temp(j) <  e_guess) THEN
           ref_walk = j
           e_guess = temp(j)
        END IF
     END DO
     temp(ref_walk) = 1.0d+30
     guess_pointers(ref_walk) = i
  END DO
!
  IF ( type_calculation == 'linear_system') THEN
       i_tri = guess_size * ( guess_size + 1 ) / 2       
       ALLOCATE( h_mat_tri_z(i_tri) )
       h_mat_tri_z = (0.d0,0.d0)
       DO i = 1, non_zero
          i_guess = guess_pointers(ibuf(1,i))
          IF (i_guess <= 0 ) cycle
          j_guess = guess_pointers(ibuf(2,i))
          IF (j_guess <= 0 ) cycle
          ii = max(i_guess,j_guess)
          jj = min(i_guess,j_guess)
          i_tri = ii * ( ii - 1 ) / 2 + jj      
          h_mat_tri_z(i_tri) = h_mat_tri_z(i_tri)  + h_buf_z(i)
       END DO
       DO  i=1,matrix_size
           i_guess = guess_pointers(i)
           IF (i_guess <= 0) cycle
           i_tri = i_guess * ( i_guess + 1 ) / 2
           h_mat_tri_z(i_tri) = diag_d(i)
       END DO
       IF( prdvd(11) ) THEN
           i_tri = 0
           title='Guess Matrix'
           write(iout,2) title
           DO i=1,guess_size
              write(iout,3) i
              write(iout,4) (h_mat_tri_z(j), j=i_tri + 1, i_tri + i )
              i_tri = i_tri + i
           END DO
       END IF
       ALLOCATE(ipivot(guess_size))
       DO  i=1,matrix_size 
           i_guess = guess_pointers(i)
           IF (i_guess <= 0) cycle
           guess_rhs(i_guess) = rhs_z(i)
       END DO
       call zhpsv('u',guess_size,1,h_mat_tri_z,ipivot,guess_rhs,guess_size,info)
       h_mat_tri_z(1:matrix_size) = (0.d0,0.d0)       
       tmp = 0.0D+00
       DO  i = 1, matrix_size
           IF (guess_pointers(i) > 0) THEN
               h_mat_tri_z(i) = guess_rhs(guess_pointers(i))
               IF (ABS(h_mat_tri_z(i)) > tmp) THEN
                   tmp = ABS(h_mat_tri_z(i))
               END IF
           END IF
       END DO
       guess_rhs(1:matrix_size) = h_mat_tri_z(1:matrix_size)
       Call iosys('open guess_rhs as new',0,0,0,'Guess_RHS.out')
       Call iosys('write real Guess_RHS to Guess_RHS',2*matrix_size,guess_rhs,0,' ')
       Call iosys('read real Guess_RHS from Guess_RHS',2*matrix_size,guess_rhs,0,' ')
       IF (prdvd(12)) THEN
           title = 'Guess Solution'
           CALL prntcmn(title,guess_rhs,guess_size,1,guess_size,1,iout,'e')
       END IF
       DEALLOCATE(temp, guess_pointers, h_mat_tri_z, ipivot)
  ELSE IF( type_calculation == 'eigenvalues') THEN
       ALLOCATE( h_mat_z(guess_size,guess_size) )
       h_mat_z = (0.d0,0.d0)
       DO i = 1, non_zero
          i_guess = guess_pointers(ibuf(1,i))
          IF (i_guess <= 0 ) cycle
          j_guess = guess_pointers(ibuf(2,i))
          IF (j_guess <= 0 ) cycle
          ii = max(i_guess,j_guess)
          jj = min(i_guess,j_guess)
          h_mat_z(ii,jj) = h_mat_z(ii,jj) + h_buf_z(i)
       END DO
       DO  i=1,matrix_size
           i_guess = guess_pointers(i)
           IF (i_guess <= 0) cycle
           h_mat_z(i_guess,i_guess) = diag_d(i)
       END DO
       DO i = 1, guess_size
          DO j = 1, i
             h_mat_z(j,i) = h_mat_z(i,j)
          END DO
       END DO
       IF (prdvd(11)) THEN
           title = 'Guess Matrix'
           CALL prntcmn(title,h_mat_z,guess_size,guess_size,guess_size,guess_size,iout,'e')
        END IF
       ALLOCATE( work(5*guess_size), rwork(5*guess_size), eigv(guess_size))
       CALL zheev('v','l',guess_size,h_mat_z,guess_size,eigv,work,5*guess_size,rwork,info)
       DO  i = 1, number_of_guess_vectors
           work(1:matrix_size) = (0.d0,0.d0)
           tmp = 0.0D+00
           DO  j = 1, matrix_size
               IF (guess_pointers(j) > 0) THEN
                   work(j) = h_mat_z(guess_pointers(j),i)
                   IF (ABS(work(j)) > tmp) THEN
                       tmp = ABS(work(j))
                   END IF
               END IF
           END DO
           guess_vectors(1:matrix_size,i) = work(1:matrix_size)
       END DO
       Call iosys('write real guess_eigenvalues to guess_vectors',guess_size,eigv,0,' ')
       Call iosys('write real guess_eigenvectors to guess_vectors',                      &
                   2*guess_size*matrix_size,h_mat_d,0,' ')
       IF (prdvd(12)) THEN
           title = 'Guess Eigenvalues'
           CALL prntfmn(title,eigv,guess_size,1,guess_size,1,iout,'e')
           title = 'Guess Eigenvectors'
           CALL prntcmn(title,h_mat_z,guess_size,guess_size,guess_size,guess_size,iout,'e')
       END IF
       DEALLOCATE(temp, guess_pointers, h_mat_z, eigv, work, rwork)
  END IF
1    FORMAT (/,t15,'Results from Guess Matrix = ', i5)
2 Format(a80)
3 Format(/,5x,'Row = ',i4)
4 Format( 15x, '(',f10.5,',',f10.5,')','(',f10.5,',',f10.5,')','(',f10.5,',',f10.5,')','(',f10.5,',',f10.5,')' )
END SUBROUTINE Guess_Solution_z
!***********************************************************************
!***********************************************************************
!deck gram_schmidt_d
  SUBROUTINE Gram_Schmidt_d(v_a,v_b,thresh,size,n_a,n_b,n_out,schmidt,prnt)
!***begin prologue     gram_schmidt
!***date written       960801  (yymmdd)
!***revision date              (yymmdd)

!***keywords           gram-schmidt, orthogonalization
!***author             barry schneider(nsf)
!***source
!***purpose            gram-schmidt orthogonalization.
!***description        a set of non-orthonormal vectors, v_b, are orthogonalized
!                      to another set of vectors, v_a, using a gram-schmidt process
!                      that checks for linear dependencies. the set of vectors,
!                      v_a, are already assumed to be internally orthonormal.

!                          v_a(n,*) = input vectors
!                          v_b(n,*) = input as non-orthogonal vectors and output
!                                     as orthonormal set.
!                          thresh  = acceptance tolerance on overlap
!                          n       = dimension of vector space
!                          nstart  = beginning vector
!                          nfin    = ending vector
!                          nout    = number vectors outputted
!                          schmdt  = perform either one or two
!                                    orthonormalizations on set
!***references
!***routines called    saxpy(clams), sdot(clams), sscal(clams)

!***end prologue       gram_schmidt
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: size
  INTEGER                                  :: n_a
  INTEGER                                  :: n_b
  INTEGER                                  :: n_out
  INTEGER                                  :: n_temp
  INTEGER                                  :: n_times
  INTEGER                                  :: upper
  INTEGER                                  :: trips
  INTEGER                                  :: i
  INTEGER                                  :: j
  REAL*8, DIMENSION(size,n_a)              :: v_a
  REAL*8, DIMENSION(size,n_b)              :: v_b
  REAL*8                                   :: thresh
  LOGICAL                                  :: schmidt
  LOGICAL                                  :: prnt
  REAL*8                                   :: overlap
  REAL*8                                   :: ddot
  REAL*8                                   :: norm
!
  n_times=1
  IF(schmidt) THEN
     n_times=2
  END IF
  IF (n_a == 0 ) THEN
      DO trips=1,n_times
         n_out = 0
         DO i=1,n_b
            norm = 1.d0 / SQRT(ddot(size,v_b(:,i),1,v_b(:,i),1))            
            DO j = 1, n_out
               overlap = ddot(size,v_b(1,j),1,v_b(1,i),1)
               CALL saxpy(size,-overlap,v_b(1,j),1,v_b(1,i),1)
            END DO
            IF(norm > thresh) THEN
               v_b(:,i) = norm * v_b(:,i)
               n_out=n_out+1
               v_b(:,n_out) = v_b(:,i)
            END IF
         END DO
      END DO
  ELSE
      upper=n_b
      DO trips=1,n_times
         n_temp=0
         DO i=1,upper
            DO j=1,n_a
               overlap = ddot(size,v_a(1,j),1,v_b(1,i),1)
               CALL saxpy(size,-overlap,v_a(1,j),1,v_b(1,i),1)
            END DO
            norm=1.d0/SQRT(ddot(size,v_b(1,i),1,v_b(1,i),1))
            IF(norm > thresh) THEN
               v_b(:,i) = norm * v_b(:,i)
               n_temp=n_temp+1
               v_b(:,n_temp) = v_b(:,i)
            END IF
         END DO
         upper=n_temp
      END DO
      DO trips=1,n_times
         n_out = 0
         DO i=1,n_temp
            norm = 1.d0 / SQRT(ddot(size,v_b(:,i),1,v_b(:,i),1))            
            DO j = 1, n_out
               overlap = ddot(size,v_b(1,j),1,v_b(1,i),1)
               CALL saxpy(size,-overlap,v_b(1,j),1,v_b(1,i),1)
            END DO
            IF(norm > thresh) THEN
               v_b(:,i) = norm * v_b(:,i)
               n_out=n_out+1
               v_b(:,n_out) = v_b(:,i)
            END IF
         END DO
      END DO      
  END IF
  IF(prnt) THEN
     WRITE(iout,1) n_a, n_b
     WRITE(iout,2) n_out
  END IF
1    FORMAT(/,1X,'schmidt orthogonalization of two sets of vectors', /,1X,  &
    'set a already assumed to be orthonormal',/,1X,  &
    'number of vectors in set a = ',i4,/,1X, 'number of vectors in set b = ',i4)
2    FORMAT(/,1X,'number of set b vectors after orthogonalization '  &
    'to set a vectors = ',i4)
END SUBROUTINE Gram_Schmidt_d
!***********************************************************************
!***********************************************************************
!deck gram_schmidt_z
  SUBROUTINE Gram_Schmidt_z(v_a,v_b,thresh,size,n_a,n_b,n_out,schmidt,prnt)
!***begin prologue     gram_schmidt
!***date written       960801  (yymmdd)
!***revision date              (yymmdd)

!***keywords           gram-schmidt, orthogonalization
!***author             barry schneider(nsf)
!***source
!***purpose            gram-schmidt orthogonalization.
!***description        a set of non-orthonormal vectors, v_b, are orthogonalized
!                      to another set of vectors, v_a, using a gram-schmidt process
!                      that checks for linear dependencies. the set of vectors,
!                      v_a, are already assumed to be internally orthonormal.

!                          v_a(n,*) = input vectors
!                          v_b(n,*) = input as non-orthogonal vectors and output
!                                    as orthonormal set.
!                          thresh  = acceptance tolerance on overlap
!                          n       = dimension of vector space
!                          nstart  = beginning vector
!                          nfin    = ending vector
!                          nout    = number vectors outputted
!                          schmidt = perform either one or two
!                                    orthonormalizations on set
!***references
!***routines called    saxpy(clams), sdot(clams), sscal(clams)

!***end prologue       gram_schmidt
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: size
  INTEGER                                  :: n_a
  INTEGER                                  :: n_b
  INTEGER                                  :: n_out
  INTEGER                                  :: n_temp
  INTEGER                                  :: n_times
  INTEGER                                  :: upper
  INTEGER                                  :: trips
  INTEGER                                  :: i
  INTEGER                                  :: j
  COMPLEX*16, DIMENSION(size,n_a)          :: v_a
  COMPLEX*16, DIMENSION(size,n_b)          :: v_b
  REAL*8                                   :: thresh
  LOGICAL                                  :: schmidt
  LOGICAL                                  :: prnt
  COMPLEX*16                               :: overlap
  REAL*8                                   :: cdotc
  REAL*8                                   :: norm
!
  n_times=1
  IF(schmidt) THEN
     n_times=2
  END IF
  IF (n_a == 0 ) THEN
      DO trips=1,n_times
         n_out = 0
         DO i=1,n_b
            norm = 1.d0 / SQRT(cdotc(size,v_b(:,i),1,v_b(:,i),1))            
            DO j = 1, n_out
               overlap = cdotc(size,v_b(1,j),1,v_b(1,i),1)
               CALL caxpy(size,-overlap,v_b(1,j),1,v_b(1,i),1)
            END DO
            IF(norm > thresh) THEN
               v_b(:,i) = norm * v_b(:,i)
               n_out=n_out+1
               v_b(:,n_out) = v_b(:,i)
            END IF
         END DO
      END DO
  ELSE
      upper=n_b
      DO trips=1,n_times
         n_temp=0
         DO i=1,upper
            DO j=1,n_a
               overlap = cdotc(size,v_a(1,j),1,v_b(1,i),1)
               CALL caxpy(size,-overlap,v_a(1,j),1,v_b(1,i),1)
            END DO
            norm=1.d0/SQRT(cdotc(size,v_b(1,i),1,v_b(1,i),1))
            IF(norm > thresh) THEN
               v_b(:,i) = norm * v_b(:,i)
               n_temp=n_temp+1
               v_b(:,n_temp) = v_b(:,i)
            END IF
         END DO
         upper=n_temp
      END DO
      DO trips=1,n_times
         n_out = 0
         DO i=1,n_temp
            norm = 1.d0 / SQRT(cdotc(size,v_b(:,i),1,v_b(:,i),1))            
            DO j = 1, n_out
               overlap = cdotc(size,v_b(1,j),1,v_b(1,i),1)
               CALL caxpy(size,-overlap,v_b(1,j),1,v_b(1,i),1)
            END DO
            IF(norm > thresh) THEN
               v_b(:,i) = norm * v_b(:,i)
               n_out=n_out+1
               v_b(:,n_out) = v_b(:,i)
            END IF
         END DO
      END DO      
  END IF
  IF(prnt) THEN
     WRITE(iout,1) n_a, n_b
     WRITE(iout,2) n_out
  END IF
1    FORMAT(/,1X,'schmidt orthogonalization of two sets of vectors', /,1X,  &
    'set a already assumed to be orthonormal',/,1X,  &
    'number of vectors in set a = ',i4,/,1X, 'number of vectors in set b = ',i4)
2    FORMAT(/,1X,'number of set b vectors after orthogonalization '  &
    'to set a vectors = ',i4)
END SUBROUTINE Gram_Schmidt_z
!***********************************************************************
!***********************************************************************
!**begin prologue     Linear_System_Driver_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           Davidson, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Davidson
!**purpose            
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Linear_System_Driver_d
!***********************************************************************
  SUBROUTINE Linear_System_Driver_d(v_in,v_out,rhs)
  IMPLICIT NONE
  REAL*8, DIMENSION(matrix_size)      :: v_in
  REAL*8, DIMENSION(matrix_size)      :: v_out
  REAL*8, DIMENSION(matrix_size)      :: rhs
  CHARACTER (LEN=5)                   :: itoc
  INTEGER                             :: i
  REAL*8                              :: ddot
!
  convergence = 'unconverged'
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
!
!
!           Initial vector
!
  local_time=secnds(0.0)
!
! normalize first vector take here as v_in
!
  anorm = 1.d0 / sqrt (ddot(matrix_size , v_in(:) , 1, v_in(:) , 1 ))
  vec_d(:,0) = anorm * v_in(:)
  total_time(1) = secnds(0.0) - local_time + total_time(1)

!
! Start iterating
!
  DO i=0, maximum_number_of_davidson_vectors
!
     WRITE(iout,3) i
!
!            Lets print the  next vector if requested
!
     IF(prdvd(2).or.dvdall) THEN
        title='Vector iteration = '//itoc(i)
        call prntfmn(title,vec_d(:,i),matrix_size,1,matrix_size,           &
                    maximum_number_of_davidson_vectors,iout,'e')
     END IF
!
!    This is the major routine which produces the solution.
!
     Call Solution_Vector(v_out,rhs,i)
!
!
!   If the number of iterations has reached the size of the
!   original matrix without breakdown, we have to be converged.
!
    IF( convergence == 'maximal') THEN
        WRITE(iout,4)
        v_in(:) = v_out(:)
        total_number_of_iterations = total_number_of_iterations + i + 1
        return
!
!  If the residual has fallen below the convergence criterion, the quit.
!
    ELSE IF(convergence == 'converged') THEN
        WRITE(iout,5) i+1  
        v_in(:) = v_out(:) 
        total_number_of_iterations = total_number_of_iterations + i + 1
        return
    END IF    
!
!   We need more to keep iterating.
!
    total_time(8) = secnds(0.0) - local_time + total_time(8)
  END DO   
  write(iout,6) maximum_number_of_davidson_vectors
  stop
1 FORMAT('***********************************************'  &
         '*************************') 
2 FORMAT(/,20X,'Beginning Davidson Iterations At Zero')
3 FORMAT(/,20x,'Davidson Iteration = ',i4)
4 FORMAT(/,20x,'Number of Iterations Maximal.  Energy = ', e15.8)
5 FORMAT(/,20x,'Convergence After ',i5, ' Iterations')
6 FORMAT(/,20x,'No Convergence After ',i5, 'Iterations')
END SUBROUTINE Linear_System_Driver_d
!***********************************************************************
!***********************************************************************
!**begin prologue     Linear_System_Driver_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           Davidson, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Davidson
!**purpose            
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Linear_System_Driver_z
!***********************************************************************
  SUBROUTINE Linear_System_Driver_z(v_in,v_out,rhs)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(matrix_size)      :: v_in
  COMPLEX*16, DIMENSION(matrix_size)      :: v_out
  COMPLEX*16, DIMENSION(matrix_size)      :: rhs
  CHARACTER (LEN=5)                       :: itoc
  INTEGER                                 :: i
  REAL*8                                  :: cdotc
!
  convergence = 'unconverged'
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
!
!
!           Initial vector
!
  local_time=secnds(0.0)
!
! normalize first vector take here as v_in
!
  anorm = 1.d0 / sqrt (cdotc(matrix_size , v_in(:) , 1, v_in(:) , 1 ))
  vec_z(:,0) = anorm * v_in(:)
  v_in(:) = anorm * vec_d(:,0)
  total_time(1) = secnds(0.0) - local_time + total_time(1)
!
! Start iterating
!
  DO i=0, maximum_number_of_davidson_vectors
!
     WRITE(iout,3) i
!
!            Lets print the  next vector if requested
!
     IF(prdvd(2).or.dvdall) THEN
        title='Vector iteration = '//itoc(i)
        call prntcmn(title,vec_z(:,i),matrix_size,1,matrix_size,           &
                    maximum_number_of_davidson_vectors,iout,'e')
     END IF
!
!    This is the major routine which produces the solution.
!
     Call Solution_Vector(v_out,rhs,i)
!
!
!   If the number of iterations has reached the size of the
!   original matrix without breakdown, we have to be converged.
!
    IF( convergence == 'maximal') THEN
        WRITE(iout,4) 
        v_in(:) =v_out(:)
        total_number_of_iterations = total_number_of_iterations + i + 1
        return
!
!           We are actually converged to the desired tolerance.
!
    ELSE IF(convergence == 'converged') THEN
        WRITE(iout,5) i+1  
        v_in(:) =v_out(:) 
        total_number_of_iterations = total_number_of_iterations + i + 1
        return
    END IF    
!
!           We need more iterations
!
    total_time(8) = secnds(0.0) - local_time + total_time(8)
  END DO   
  write(iout,6) maximum_number_of_davidson_vectors
  stop
1 FORMAT('***********************************************'  &
         '*************************') 
2 FORMAT(/,20X,'Beginning Davidson Iterations At Zero')
3 FORMAT(/,20x,'Davidson Iteration = ',i4)
4 FORMAT(/,20x,'Number of Iterations Maximal.  Energy = ', e15.8)
5 FORMAT(/,20x,'Convergence After ',i5, ' Iterations')
6 FORMAT(/,20x,'No Convergence After ',i5, 'Iterations')
END SUBROUTINE Linear_System_Driver_z
!***********************************************************************
!***********************************************************************
!**begin prologue     Solution_Vector_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           Davidson
!**
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            Perform a Davidson iteration.
!**                   Here we solve the projected linear system,                   
!***                             A|X> = |B> 
!***                  Step one: Multiply the input vector by A
!***                  Step two: Form the projected linear system and solve.
!***                  Step three: Check convergence.
!***                  Step four: If converged, quit.  If not add residual to
!***                             vector space and continue.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Solution_Vector_d
!***********************************************************************
  SUBROUTINE Solution_Vector_d(v_out,rhs,it)
  IMPLICIT NONE
  REAL*8, DIMENSION(matrix_size)                                        :: v_out
  REAL*8, DIMENSION(matrix_size)                                        :: rhs
  INTEGER                                                               :: i
  INTEGER                                                               :: j
  INTEGER                                                               :: it
  INTEGER                                                               :: ii
  INTEGER                                                               :: index
  INTEGER                                                               :: info
  REAL*8                                                                :: ddot
  CHARACTER(LEN=4)                                                      :: itoc
!
!            Note that it runs from 0 to maxvec
!            The starting vector is assumed to be normalized.
!
!            Lets compute the value of A on the input vector in order to
!            begin the generation of the next vector.
!
  local_time=secnds(0.0)
  call packed_symmetric_matrix_on_vector                                  &
                                       (h_buf_d,                          &
                                        ibuf,                             &
                                        diag_d,                           &
                                        vec_d(1:matrix_size,it),          &
                                        h_vectors_d(1:matrix_size,it),    &
                                        non_zero)
  IF(prdvd(3).or.dvdall) THEN
     title='H_on_Vector iteration = '//itoc(it)
     call prntfmn(title,h_vectors_d(:,it),matrix_size,1,matrix_size,1,iout,'e')
  END IF
  total_time(2) = secnds(0.0) - local_time + total_time(2)
!
!     Calculate the new elements, right hand side and solve small linear system
! 
!
  Call Solve_Small_Linear_System(h_mat_tri_d,h_mat_work_tri_d,small_rhs_d,small_rhs_work_d,it)
!
  total_time(3) = secnds(0.0) - local_time + total_time(3)
!
!             Test convergence
!
  local_time=secnds(0.0)
  call Convergence_Test(v_out,vec_d,h_vectors_d,small_rhs_work_d,rhs,it)
  total_time(4) = secnds(0.0) - local_time + total_time(4)
!
1 Format(a80)
2 Format(/,5x,'Row = ',i4)
3 Format( (15x,5f10.5) )
END SUBROUTINE Solution_Vector_d
!***********************************************************************
!***********************************************************************
!**begin prologue     Solution_Vector_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           Davidson
!**
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            Perform a Davidson iteration.
!**                   Here we solve the projected linear system,                   
!***                             A|X> = |B> 
!***                  Step one: Multiply the input vector by A
!***                  Step two: Form the projected linear system and solve.
!***                  Step three: Check convergence.
!***                  Step four: If converged, quit.  If not add residual to
!***                             vector space and continue.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Solution_Vector_z
!***********************************************************************
  SUBROUTINE Solution_Vector_z(v_out,rhs,it)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(matrix_size)                                    :: v_out
  COMPLEX*16, DIMENSION(matrix_size)                                    :: rhs
  INTEGER                                                               :: i
  INTEGER                                                               :: j
  INTEGER                                                               :: it
  INTEGER                                                               :: ii
  INTEGER                                                               :: index
  INTEGER                                                               :: info
  COMPLEX*16                                                            :: cdotc
  CHARACTER(LEN=4)                                                      :: itoc
!
!            Note that it runs from 0 to maxvec
!            The starting vector is assumed to be normalized.
!
!            Lets compute the value of A on the input vector in order to
!            begin the generation of the next vector.
!
  local_time=secnds(0.0)
  call packed_symmetric_matrix_on_vector                                  &
                                       (h_buf_z,                          &
                                        ibuf,                             &
                                        diag_z,                           &
                                        vec_z(1:matrix_size,it),          &
                                        h_vectors_z(1:matrix_size,it),    &
                                        non_zero)
  IF(prdvd(3).or.dvdall) THEN
     title='H_on_Vector iteration = '//itoc(it)
     call prntcmn(title,h_vectors_z(:,it),matrix_size,1,matrix_size,1,iout,'e')
  END IF
  total_time(2) = secnds(0.0) - local_time + total_time(2)
!
!     Calculate the new elements, right hand side and solve small linear system
! 

!
!             Compute solution to the projected linear system
!
  Call Solve_Small_Linear_System(h_mat_tri_d,h_mat_work_tri_d,small_rhs_d,small_rhs_work_d,it)
!
  total_time(3) = secnds(0.0) - local_time + total_time(3)
!
!             Test convergence
!
  local_time=secnds(0.0)
  call Convergence_Test(v_out,vec_z,h_vectors_z,small_rhs_work_z,rhs,it)
  total_time(4) = secnds(0.0) - local_time + total_time(4)
!
1 Format(a80)
2 Format(/,5x,'Row = ',i4)
3 Format( 15x, '(',f10.5,',',f10.5,')','(',f10.5,',',f10.5,')',             &
               '(',f10.5,',',f10.5,')','(',f10.5,',',f10.5,')' )
END SUBROUTINE Solution_Vector_z
!***********************************************************************
!***********************************************************************
!deck Solve_Small_Linear_System_d
!***begin prologue     Solve_Small_Linear_System_d
!***date written       980420   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Solve Small davidson matrix
!***author             schneider, barry (nsf)
!***source
!***
!***references
!***routines called
!***end prologue       Solve_Small_Linear_System_d
  SUBROUTINE Solve_Small_Linear_System_d(h_mat,h_mat_work,small_rhs,small_rhs_work,it)
  REAL*8, DIMENSION(0:n_tri)                                            :: h_mat
  REAL*8, DIMENSION(0:n_tri)                                            :: h_mat_work
  REAL*8, DIMENSION(0:maximum_number_of_davidson_vectors)               :: small_rhs
  REAL*8, DIMENSION(0:maximum_number_of_davidson_vectors)               :: small_rhs_work
  REAL*8                                                                :: ddot
  INTEGER                                                               :: i
  INTEGER                                                               :: j
  INTEGER                                                               :: ii
  INTEGER                                                               :: it
  INTEGER                                                               :: info
  INTEGER                                                               :: index
  CHARACTER(LEN=4)                                                      :: itoc
!
!
!     Calculate the new elements and right hand side
! 
  ii = it * ( it + 1 ) /2
  DO i = 0, it
     index = ii + i
     h_mat(index) = ddot(matrix_size,vec_d(:,i),1,h_vectors_d(:,it),1)
  END DO
  h_mat_work(0:index) = h_mat(0:index)
  small_rhs(it) = ddot(matrix_size,vec_d(:,it),1,rhs_d(:),1)
  small_rhs_work(0:it) = small_rhs(0:it)
  IF(prdvd(11).or.dvdall) THEN
     title='Small Matrix iteration = '//itoc(it)
     write(iout,1) title
     ii = 0
     DO i=0, it
        write(iout,2) i
        write(iout,3) (h_mat_work(j), j=ii, ii + i )
        ii = ii + i + 1
     END DO
     title='Small RHS iteration = '//itoc(it)
     call prntfmn(title,small_rhs_work,it+1,1,it+1,1,iout,'e')
  END IF
!
  CALL dspsv('u',it+1,1,h_mat_work,ipvt,small_rhs_work,it+1,info)
  IF(info /= 0) THEN
     Write(iout,*) '     Warning Singular System at This Step'
  END IF
  IF(prdvd(11).or.dvdall) THEN
     title='Small Matrix Solution iteration = '//itoc(it)
     call prntfmn(title,small_rhs_work,it+1,1,it+1,1,iout,'e')
  END IF
!
1 Format(a80)
2 Format(/,5x,'Row = ',i4)
3 Format( (15x,5f10.5) )
END SUBROUTINE Solve_Small_Linear_System_d
!***********************************************************************
!***********************************************************************
!deck Solve_Small_Linear_System_z
!***begin prologue     Solve_Small_Linear_System_z
!***date written       980420   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           Solve Small davidson matrix
!***author             schneider, barry (nsf)
!***source
!***
!***references
!***routines called
!***end prologue       Solve_Small_Linear_System_z
  SUBROUTINE Solve_Small_Linear_System_z(h_mat,small_rhs,it)
  COMPLEX*16, DIMENSION(0:n_tri)                                        :: h_mat
  COMPLEX*16, DIMENSION(0:n_tri)                                        :: h_mat_work
  COMPLEX*16, DIMENSION(0:maximum_number_of_davidson_vectors)           :: small_rhs
  COMPLEX*16, DIMENSION(0:maximum_number_of_davidson_vectors)           :: small_rhs_work
  COMPLEX*16                                                            :: cdotc
  INTEGER                                                               :: it  
  INTEGER                                                               :: i
  INTEGER                                                               :: j
  INTEGER                                                               :: index
  INTEGER                                                               :: ii
  INTEGER                                                               :: info
  CHARACTER(LEN=4)                                                      :: itoc
!
  ii = it * ( it + 1 ) /2
  DO i = 0, it
     index = ii + i
     h_mat(index) = cdotc(matrix_size,vec_z(:,i),1,h_vectors_z(:,it),1)
  END DO
  h_mat_work(0:index) = h_mat(0:index)
  small_rhs(it) = cdotc(matrix_size,vec_z(:,it),1,rhs_z(:),1)
  small_rhs_work(0:it) = small_rhs(0:it)
  IF(prdvd(11).or.dvdall) THEN
     title='Small Matrix iteration = '//itoc(it)
     ii = 0
     DO i=0, it
        write(iout,2) i
        write(iout,3) (h_mat_work(j), j=ii, ii + i )
        ii = ii + i + 1
     END DO
     title='Small RHS iteration = '//itoc(it)
     call prntcmn(title,small_rhs_work,it+1,1,it+1,1,iout,'e')
  END IF
  CALL zhpsv('u',it+1,1,h_mat_work,ipvt,small_rhs_work,it+1,info)
  IF(info /= 0) THEN
     Write(iout,*) '     Warning Singular System at This Step'
  END IF
  IF(prdvd(11).or.dvdall) THEN
     title='Small Matrix Solution iteration = '//itoc(it)
     call prntcmn(title,small_rhs_work,it+1,1,it+1,1,iout,'e')
  END IF
!
1 Format(a80)
2 Format(/,5x,'Row = ',i4)
3 Format( 15x, '(',f10.5,',',f10.5,')','(',f10.5,',',f10.5,')','(',f10.5,',',f10.5,')','(',f10.5,',',f10.5,')' )
END SUBROUTINE Solve_Small_Linear_System_z
!***********************************************************************
!***********************************************************************
!deck Solve_Small_Eigenvalue_Problem_d
!***begin prologue     Solve_Small_Eigenvalue_Problem_d
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           small davidson matrix
!***author             schneider, barry (nsf)
!***source
!***
!***references
!***routines called
!***end prologue       Solve_Small_Eigenvalue_Problem_d
  SUBROUTINE Solve_Small_Eigenvalue_Problem_d(h_mat,h_mat_work,eigen_val_work,it,  &  
                                              first_vector,last_vector)
  IMPLICIT NONE
  REAL*8, DIMENSION(0:maximum_number_of_davidson_vectors,0:last_vector)            &
                                                               :: h_mat
  REAL*8, DIMENSION(0:maximum_number_of_davidson_vectors,0:last_vector)            &
                                                               :: h_mat_work
  REAL*8, DIMENSION(0:maximum_number_of_davidson_vectors)      :: eigen_val_work
  INTEGER                                                      :: first_vector
  INTEGER                                                      :: last_vector
  INTEGER                                                      :: it
  INTEGER                                                      :: i
  INTEGER                                                      :: j
  INTEGER                                                      :: ii
  INTEGER                                                      :: index
  INTEGER                                                      :: info
  REAL*8                                                       :: ddot
  CHARACTER (LEN=3)                                            :: itoc
  DO i=first_vector, last_vector
     DO j = 0, i
        h_mat(i,j) = ddot(matrix_size,vec_d(:,i),1,h_vectors_d(:,j),1)
        h_mat(j,i) = h_mat(i,j)
     END DO
  END DO
  h_mat_work(0:last_vector,0:last_vector) =  h_mat(0:last_vector,0:last_vector)
  IF(prdvd(11).or.dvdall) THEN
     title='Small Matrix Iteration = '//itoc(it)
     Call prntfm(title,h_mat,last_vector+1,last_vector+1,                          &
                 maximum_number_of_davidson_vectors+1,                             &
                 maximum_number_of_davidson_vectors+1,iout)
  END IF
  CALL dsyev('v','l',last_vector+1,h_mat_work,                                     &
              maximum_number_of_davidson_vectors+1,eigen_val_work,work_d,          &
              5*(maximum_number_of_davidson_vectors+1),info)
  IF(prdvd(11).or.dvdall) THEN
     title='Eigenvalues = '//itoc(it)
     Call prntfm(title,eigen_val_work,last_vector+1,1,last_vector+1,1,iout)
  END IF
  END SUBROUTINE Solve_Small_Eigenvalue_Problem_d
!***********************************************************************
!***********************************************************************
!deck Solve_Small_Eigenvalue_Problem_z
!***begin prologue     Solve_Small_Eigenvalue_Problem_z
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           small davidson matrix
!***author             schneider, barry (nsf)
!***source
!***
!***references
!***routines called
!***end prologue       Solve_Small_Eigenvalue_Problem_z
  SUBROUTINE Solve_Small_Eigenvalue_Problem_z(h_mat,h_mat_work,eigen_val_work,it,     &
                                              first_vector,last_vector)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(0:maximum_number_of_davidson_vectors,0:last_vector)           &
                                                               :: h_mat
  COMPLEX*16, DIMENSION(0:maximum_number_of_davidson_vectors,0:last_vector)           &
                                                               :: h_mat_work
  REAL*8, DIMENSION(0:maximum_number_of_davidson_vectors)      :: eigen_val_work
  INTEGER                                                      :: first_vector
  INTEGER                                                      :: last_vector
  INTEGER                                                      :: it
  INTEGER                                                      :: i
  INTEGER                                                      :: j
  INTEGER                                                      :: ii
  INTEGER                                                      :: index
  INTEGER                                                      :: info
  COMPLEX*16                                                   :: cdotc
  CHARACTER (LEN=3)                                            :: itoc
  DO i=first_vector, last_vector
     DO j = 0, i
        h_mat(i,j) = cdotc(matrix_size,vec_z(:,i),1,h_vectors_z(:,j),1)
        h_mat(j,i) = conjg( h_mat(i,j) )
     END DO
  END DO
  h_mat_work(0:last_vector,0:last_vector) = h_mat(0:last_vector,0:last_vector)
  IF(prdvd(11).or.dvdall) THEN
     title='Small Matrix Iteration = '//itoc(it)
     Call prntcm(title,h_mat,last_vector+1,last_vector+1,                              &
                 maximum_number_of_davidson_vectors+1,                                 &
                 maximum_number_of_davidson_vectors+1,iout)
  END IF
  CALL zheev('v','l',last_vector+1,h_mat_work,                                         &
                      maximum_number_of_davidson_vectors+1,eigen_val_work,             &
                      work_z,5*(maximum_number_of_davidson_vectors+1),                 &
                      rwork,info)
  IF(prdvd(11).or.dvdall) THEN
     title='Eigenvalues = '//itoc(it)
     Call prntfm(title,eigen_val_work,last_vector+1,1,last_vector+1,1,iout)
  END IF
  END SUBROUTINE Solve_Small_Eigenvalue_Problem_z
!***********************************************************************
!***********************************************************************
!**begin prologue     Convergence_Test_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Davidson
!**purpose            
!***                  
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Convergence_Test_d
!***********************************************************************
  SUBROUTINE Convergence_Test_d(v_out,vectors,h_vectors,small_rhs,rhs,it)
  IMPLICIT NONE
  REAL*8, DIMENSION(matrix_size)                                        :: v_out
  REAL*8, DIMENSION(matrix_size,0:maximum_number_of_davidson_vectors)   :: vectors
  REAL*8, DIMENSION(matrix_size,0:maximum_number_of_davidson_vectors)   :: h_vectors
  REAL*8, DIMENSION(0:maximum_number_of_davidson_vectors)               :: small_rhs
  REAL*8, DIMENSION(matrix_size)                                        :: rhs
  INTEGER                                                               :: it
  INTEGER                                                               :: i
  INTEGER                                                               :: n_out
  REAL*8                                                                :: ddot
  REAL*8                                                                :: ERR
  REAL*8                                                                :: test
!
!
!     Calculate the residual rhs = rhs - AX
!
  v_out(:) = rhs(:)
  Call ambcxx(v_out,h_vectors,small_rhs,matrix_size,it+1,1,matrix_size,  &
              matrix_size,maximum_number_of_davidson_vectors)
!
!     Test Convergence
!
  ERR = SQRT (ddot(matrix_size,v_out,1,v_out,1) )
!
 IF(prdvd(14).or.dvdall) THEN
     title='residual'
     call prntfmn(title,v_out,matrix_size,1,matrix_size,1,iout,'e')
 END IF
!
!
  IF ( ERR > davidson_convergence.and.it == maximum_number_of_davidson_vectors) THEN
      convergence = 'maximal'
      return
  ELSE IF (ERR <= davidson_convergence) THEN
    convergence='converged'
    Call ebcxx(v_out,vectors,small_rhs,matrix_size,it+1,1,matrix_size,  &
               matrix_size,maximum_number_of_davidson_vectors)
    title='Final Solution'
    call prntfmn(title,v_out,matrix_size,1,matrix_size,1,iout,'e')
  ELSE
    convergence='unconverged'
    IF (preconditioner == 'diagonal' ) THEN
        DO i = 1, matrix_size
           test =ABS(diag_d(i))
           IF(ABS(test) >= near_zero) THEN
               v_out(i) = v_out(i)/test
           ELSE
               v_out(i) = one
           END IF
        END DO
    ELSE IF(preconditioner == 'other' ) THEN
        Call lnkerr('not yet implimented')
    END IF
    CALL gram_schmidt_d(vectors,v_out,overlap_tolerance,matrix_size,it+1,1,n_out,.true.,.false.)
    IF (n_out /= 1) THEN
        write(iout,*) 'No More Vectors.  Quit'
        stop
    END IF
    vectors(:,it+1) = v_out(:)
    IF(prdvd(1).or.dvdall) THEN
       title='basis'
       call prntfmn(title,vectors,matrix_size,it+2,matrix_size,maximum_number_of_davidson_vectors+1,iout,'e')
    END IF
  END IF
  WRITE(iout,1) ERR, convergence
1 FORMAT(/,10x,'RMS Error       = ',e15.8,/,10X, 'Status          = ',a16)
END SUBROUTINE Convergence_Test_d
!***********************************************************************
!***********************************************************************
!**begin prologue     Convergence_Test_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Davidson
!**purpose            
!***                  
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Convergence_Test_z
!***********************************************************************
  SUBROUTINE Convergence_Test_z(v_out,vectors,h_vectors,small_rhs,rhs,it)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(matrix_size)                                        :: v_out
  COMPLEX*16, DIMENSION(matrix_size,0:maximum_number_of_davidson_vectors)   :: vectors
  COMPLEX*16, DIMENSION(matrix_size,0:maximum_number_of_davidson_vectors)   :: h_vectors
  COMPLEX*16, DIMENSION(0:maximum_number_of_davidson_vectors)               :: small_rhs
  COMPLEX*16, DIMENSION(matrix_size)                                        :: rhs
  INTEGER                                                                   :: it
  INTEGER                                                                   :: i
  INTEGER                                                                   :: n_out
  COMPLEX*16                                                                :: cdotc
  REAL*8                                                                    :: ERR
  REAL*8                                                                    :: test
!
!
!     Calculate the residual rhs = rhs - AX
!
  v_out(:) = rhs(:)
  Call cambcx(v_out,matrix_size,h_vectors,matrix_size,small_rhs,maximum_number_of_davidson_vectors,  &
              matrix_size,it+1,1)
  ERR = SQRT (cdotc(matrix_size,v_out,1,v_out,1) )
!
  WRITE(iout,1) ERR, convergence
  IF ( ERR > davidson_convergence.and.it == maximum_number_of_davidson_vectors) THEN
      convergence = 'maximal'
      return
  ELSE IF(ERR <= davidson_convergence) THEN
    convergence='converged'
    Call cebcx(v_out,matrix_size,vectors,matrix_size,small_rhs,maximum_number_of_davidson_vectors,matrix_size,it+1,1)
    title='Final Solution'
    call prntcmn(title,v_out,matrix_size,1,matrix_size,1,iout,'e')
  ELSE
    convergence='unconverged'
    IF (preconditioner == 'diagonal' ) THEN
        DO i = 1, matrix_size
           test =ABS(diag_z(i))
           IF(ABS(test) >= near_zero) THEN
               v_out(i) = v_out(i)/test
           ELSE
               v_out(i) = one
           END IF
        END DO
    ELSE IF(preconditioner == 'other' ) THEN
        Call lnkerr('not yet implimented')
    END IF
    CALL gram_schmidt_z(vectors,v_out,overlap_tolerance,matrix_size,it+1,1,n_out,.true.,.false.)
    IF (n_out /= 1) THEN
        write(iout,*) 'No More Vectors.  Quit'
        stop
    END IF
    vectors(:,it+1) = v_out(:)
    IF(prdvd(1).or.dvdall) THEN
       title='basis'
       call prntcmn(title,vectors,matrix_size,it+2,matrix_size,maximum_number_of_davidson_vectors+1,iout,'e')
    END IF
  END IF
1 FORMAT(/,10x,'RMS Error       = ',e15.8,/,5X, 'Status          = ',a16)
END SUBROUTINE Convergence_Test_z
!***********************************************************************
!***********************************************************************
!**begin prologue     Initialize_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Davidson
!**purpose            
!***                  
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Initialize_d
!***********************************************************************
  SUBROUTINE Initialize_d(vectors_out,scr,number_in,number_out,number_converged)
  IMPLICIT NONE
  REAL*8, DIMENSION(matrix_size,number_out)         :: vectors_out
  REAL*8, DIMENSION(matrix_size,number_out)         :: scr
  INTEGER                                           :: number_in
  INTEGER                                           :: number_out
  INTEGER                                           :: number_converged
  INTEGER                                           :: actual_number
!
!
  IF(number_converged /= 0) THEN
     Call iosys('read real converged_davidson_vectors from davidson_vectors', &
                 number_converged*matrix_size,scr,0, ' ')
     CALL gram_schmidt(scr,vectors_out,thresh,matrix_size,                      &
                       number_converged,number_in,actual_number,.true.,.false.)
  END IF
  IF(actual_number /= 0) THEN
     CALL gram_schmidt(scr,vectors_out,thresh,matrix_size,0,actual_number,      &
                       number_out,.true.,.false.)
  END IF
!
! If we have no more vectors, we have to quit.
!
  IF(number_out == 0) THEN
     WRITE(iout,1)
     CALL lnkerr('quit davidson. no more trial vectors '// 'possible')
  END IF
1 FORMAT(/,5X,'cannot even begin davidson calculation:',/,5X,  &
              'orthonormalization of initial vectors yields null ' 'set')
  END SUBROUTINE Initialize_d
!***********************************************************************
!***********************************************************************
!**begin prologue     Initialize_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Davidson
!**purpose            
!***                  
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Initialize_z
!***********************************************************************
  SUBROUTINE Initialize_z(vectors_out,scr,number_in,number_out,number_converged)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(matrix_size,number_out)     :: vectors_out
  COMPLEX*16, DIMENSION(matrix_size,number_out)     :: scr
  INTEGER                                           :: number_in
  INTEGER                                           :: number_out
  INTEGER                                           :: number_converged
  INTEGER                                           :: actual_number
!
!
  IF(number_converged /= 0) THEN
     Call iosys('read real converged_davidson_vectors from davidson_vectors', &
                 2*number_converged*matrix_size,scr,0, ' ')
     CALL gram_schmidt(scr,vectors_out,thresh,matrix_size,                      &
                       number_converged,number_in,actual_number,.true.,.false.)
  END IF
  IF(actual_number /= 0) THEN
     CALL gram_schmidt(scr,vectors_out,thresh,matrix_size,0,actual_number,      &
                       number_out,.true.,.false.)
  END IF
!
! If we have no more vectors, we have to quit.
!
  IF(number_out == 0) THEN
     WRITE(iout,1)
     CALL lnkerr('quit davidson. no more trial vectors '// 'possible')
  END IF
1 FORMAT(/,5X,'cannot even begin davidson calculation:',/,5X,  &
              'orthonormalization of initial vectors yields null ' 'set')
  END SUBROUTINE Initialize_z
!***********************************************************************
!***********************************************************************
           END MODULE Davidson_Module
!***********************************************************************
!***********************************************************************



