!***********************************************************************
! Generalized_Lanczos_Module
!**begin prologue     Generalized_Lanczos_Module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains all of the major subroutines to propagate
!***                  a wavefunction in time using the Lanczos
!***                  algorithm.  The routine works for standard and generalized 
!***                  symmetric problems.
!***                  Explicit interfaces are used to allow a transparent use of 
!***                  generic subroutines which work for both real and complex vectors.  
!***                  This feature permits a single code to be used for both real and
!***                  imaginary time propagation.
!***description       Given a starting vector, a number of iterations
!***                  are performed until the time propagated solution
!***                  satisfies a fixed accuracy criterion.  The criterion
!***                  used depends on whether one is propagating in real or
!***                  imaginary time.  In imaginary time, the number of 
!***                  iterations at a given step is controlled by the 
!***                  convergence of the eigenvalue of interest.  The
!***                  converged vector provides the starting vector for the
!***                  next time step.  For real time, the number of iterations
!***                  at a given step depends on the RMS deviation of the 
!***                  wavefunction.  When that deviation is smaller than some
!***                  preset value, the iterations are terminated.  The 
!***                  converged wavefunction is used as the starting wavefunction
!***                  at the next time step.
!***       
!***                                  Details
!***
!***                  Perform the Lanczos recursion at a given timestep.
!***                  b_{i+1} S |V_{i+1}> = [ H - a_{i} S ] |V_{i}> 
!                                               -   b_{i} S |V_{i-1}>
!                                        T
!                                   S = U U
!                                    -1     -1  -T
!                                   S  = ( U   U   )
!***                  We rewite this as,
!***                                         -T  -1
!***                  b_{i+1} |X_{i+1}> = [ U H U - a_{i} ] |X_{i}> 
!***                                               -   b_{i} |X_{i-1}>
!***                          
!***                                |X_{i}> = U |V_{i}>
!***
!***                  The Lanczos recursion for a generalized eigenvalue problem requires
!***                  the factoization of the S matrix and the solution of two triangular
!***                  linear systems at each iteration ans well as a matrix multiply.
!***                  The multiply and linear system solves are all comparable in flop count.
!***references
!***modules needed    See USE statements below
!***comments          In this portable version I have disabled all unnecessary
!***                  writing to files.  The original Fortran is commented out.
!***                  In addition, there is no option to compute the autocorrelation
!***                  function as this would require reading and manipulating the
!***                  initial state wavefunction from a file.
!***end prologue      Generalized_Lanczos_Module
!***********************************************************************
!***********************************************************************
                           MODULE Generalized_Lanczos_Module
                           USE Global_Time_Propagation_Module, ONLY:packed_matrices, matrix_type, non_orth
                           USE full_matrix_vector_multiply_module
                           USE packed_matrix_vector_multiply_module
                           USE dvr_matrix_vector_multiply_module
                           USE Preconditioner_Module
                           USE Lanczos_Global,  ONLY: lanczos_convergence,              &
                                                maximum_number_of_time_subintervals,    &
                                                total_number_of_iterations,             &
                                                t_start,                                &
                                                t_end,                                  &
                                                save_deltat
                           USE Exponential_on_Vector_Module
                           IMPLICIT NONE
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                           INTERFACE Generalized_Lanczos_Matmul
             MODULE PROCEDURE Generalized_Lanczos_Matmul_d,                             &
                              Generalized_Lanczos_Matmul_z
                       END INTERFACE Generalized_Lanczos_Matmul
!
                           INTERFACE Generalized_Lanczos
             MODULE PROCEDURE Generalized_Lanczos_d,                                    &
                              Generalized_Lanczos_z
                       END INTERFACE Generalized_Lanczos
!
                           INTERFACE Generalized_Lanczos_Vec
             MODULE PROCEDURE Generalized_Lanczos_Vec_d,                                &
                              Generalized_Lanczos_Vec_z
                       END INTERFACE Generalized_Lanczos_Vec
!
                           INTERFACE Generalized_Lanczos_Hamiltonian
             MODULE PROCEDURE Generalized_Lanczos_Hamiltonian
                       END INTERFACE Generalized_Lanczos_Hamiltonian
!
                           INTERFACE Generalized_Lanczos_Convergence_Test
             MODULE PROCEDURE Generalized_Lanczos_convergence_test_d,                   &
                              Generalized_Lanczos_Convergence_Test_z  
                       END INTERFACE Generalized_Lanczos_Convergence_Test
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!**begin prologue     Generalized_Lanczos_Matmul_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Perform the matrix vector multiply for both the standard
!***                  generalized problem.
!**description        
!**                   
!**references
!**routines called    
!**end prologue       Generalized_Lanczos_Matmul_d
!***********************************************************************
  SUBROUTINE Generalized_Lanczos_Matmul_d(v_i,v_o,it)
  IMPLICIT NONE
  TYPE(REAL_PROP_VAR)                               :: prop_mat_d
  TYPE(REAL_MATRICES)                               :: mat_d
  REAL(idp), INTENT(in),  DIMENSION(:)              :: v_i
  REAL(idp), INTENT(out), DIMENSION(:)              :: v_o
  REAL(idp)                                         :: ddot
  INTEGER                                           :: it 
!
!           In this section of code we compute;
!
!           -T   -1
!         U   H U  |v_i> = M |v_i> 
!
!           If the problem is not a generalized problem, its easier but the coding here
!           allows for both cases as well as packed and non-packed matrices.
!
!  write(iout,*) v_i 
  IF (non_orth) THEN  ! This is for the real generalized problem when the Hamiltonian is
                      ! is not transformed to a standard problem
      IF (packed_matrices) THEN ! The routines called are specialized for matrices where the
                                ! elements have been packed without any zeros.  The actual routines
                                ! called are in a module that is in the library.
!
!                             -1      
!            Solve  |v_a> =  U  |v_i>

          Call Packed_Solve (                                                             &
                             prop_mat_d%v_scr,                                            &
                             v_i,                                                         &
                             mat_var(1)%non_zero_columns,                                 &
                             mat_var(1)%row_index,                                        &
                             mat_var(1)%mat_d%packed_columns,                             &
                             mat_var(1)%mat_d%matrix_diagonal,                            &
                             mat_var(1)%number,'n',it)
!          write (iout,*) mat_var(1)%non_zero_columns
!          write (iout,*) mat_var(1)%row_index
!          write (iout,*) mat_var(1)%packed_columns_d
!          write (iout,*) mat_var(1)%matrix_diagonal_d
!
!            Now multiply by the Hamiltonian
!
!            |v_b> = H |v_a>
!
          Call Column_Packed_Symmetric_Matrix_on_Vector (                                 &
                                                         prop_mat_d%v_scr,                &
                                                         prop_mat_d%local_scratch,        &
                                                         mat_var(3)%mat_d%packed_columns, &
                                                         mat_var(3)%non_zero_columns,     &
                                                         mat_var(3)%row_index,            &
                                                         mat_var(3)%mat_d%matrix_diagonal )
!          write (iout,*) mat_var(3)%non_zero_columns
!          write (iout,*) mat_var(3)%row_index
!          write (iout,*) mat_var(3)%packed_columns_d
!          write (iout,*) mat_var(3)%matrix_diagonal_d
!
!                                 -T      
!            Solve  |v_o> =  U  |v_o>
!
          Call Packed_Solve (                                                             &
                             v_o,                                                         &
                             prop_mat_d%local_scratch,                                    &
                             mat_var(2)%non_zero_columns,                                 &
                             mat_var(2)%row_index,                                        &
                             mat_var(2)%mat_d%packed_columns,                             &
                             mat_var(1)%mat_d%matrix_diagonal,                            &
                             mat_var(2)%number,'t',it)
!          write (iout,*) mat_var(2)%non_zero_columns
!          write (iout,*) mat_var(2)%row_index
!          write (iout,*) mat_var(2)%packed_columns_d
!          write (iout,*) mat_var(1)%matrix_diagonal_d
!          write(iout,*) v_o
      ELSE  ! Here the matrices are not packed but the procedure is the same.
          Call Triangular_Solve(prop_mat_d%upper,                                         &
                                prop_mat_d%v_scr,                                         &
                                v_i,                                                      &
                                'n',it)
          Call Matrix_Vector_Multiply( prop_mat_d%triangle_hamiltonian,                   &
                                       prop_mat_d%v_scr,                                  &
                                       prop_mat_d%local_scratch)
          Call Triangular_Solve(prop_mat_d%upper,                                         &
                                prop_mat_d%local_scratch,                                 &
                                v_o,                                                      &
                                't',it)
      END IF
  ELSE  ! This is for the standard case either with or without packed matrices
      IF (packed_matrices) THEN 
          Call Column_Packed_Symmetric_Matrix_on_Vector (                                 &
                                                         v_i,                             &
                                                         v_o,                             &
                                                         mat_var(3)%mat_d%packed_columns, &
                                                         mat_var(3)%non_zero_columns,     &
                                                         mat_var(3)%row_index,            &
                                                         mat_var(3)%mat_d%matrix_diagonal )
      END IF
      IF (matrix_type == 'triangle') THEN
          Call Matrix_Vector_Multiply( prop_mat_d%triangle_hamiltonian,                   &
                                       v_i,                                               &
                                       v_o)
      ELSE IF (matrix_type == 'fedvr') THEN
          Call finite_element_m_v(v_i,v_o)

      END IF
  END IF
END SUBROUTINE Generalized_Lanczos_Matmul_d
!***********************************************************************
!***********************************************************************
!**begin prologue     Generalized_Lanczos_Matmul_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Same as real toutine but for complex vectors
!**description        
!**                   
!**references
!**routines called    
!**end prologue       Generalized_Lanczos_Matmul_z
!***********************************************************************
  SUBROUTINE Generalized_Lanczos_Matmul_z(v_i,v_o,it)
  IMPLICIT NONE
  TYPE(COMPLEX_PROP_VAR)                           :: prop_mat_z
  TYPE(REAL_PROP_VAR)                              :: prop_mat_d
  TYPE(REAL_MATRICES)                              :: mat_d
  COMPLEX(idp), DIMENSION(:)                       :: v_i
  COMPLEX(idp), DIMENSION(:)                       :: v_o
  INTEGER                                          :: it 
!
!           In this section of code we compute;
!
!           -T   -1
!         U   H U  |v_i>
!
!           If the problem is not a generalized problem, its easier but the coding here
!           allows for both cases as well as packed and non-packed matrices.
!
  IF (non_orth) THEN
      IF (packed_matrices) THEN 
!
!                               -1      
!            Solve  |v_o> =  U  |v_i>
!
          Call Packed_Solve (                                                             &
                             prop_mat_z%v_scr,                                            &
                             v_i,                                                         &
                             mat_var(1)%non_zero_columns,                                 &
                             mat_var(1)%row_index,                                        &
                             mat_var(1)%mat_d%packed_columns,                             &
                             mat_var(1)%mat_d%matrix_diagonal,                            &
                             mat_var(1)%number,'n',it)
!
!            Now multiply by the Hamiltonian
!
!            |v_b> - H |v_a>
!
          Call Column_Packed_Symmetric_Matrix_on_Vector (                                 &
                                                         prop_mat_z%v_scr,                &
                                                         prop_mat_z%local_scratch,        &
                                                         mat_var(3)%mat_d%packed_columns, &
                                                         mat_var(3)%non_zero_columns,     &
                                                         mat_var(3)%row_index,            &
                                                         mat_var(3)%mat_d%matrix_diagonal )
!                             -T      
!            Solve  |v_o> =  U v_b>
!
          Call Packed_Solve (                                                             &
                             v_o,                                                         &
                             prop_mat_z%local_scratch,                                    &
                             mat_var(2)%non_zero_columns,                                 &
                             mat_var(2)%row_index,                                        &
                             mat_var(2)%mat_d%packed_columns,                             &
                             mat_var(1)%mat_d%matrix_diagonal,                            &
                             mat_var(2)%number,'t',it)
      ELSE
          Call Triangular_Solve(prop_mat_d%upper,                                         &
                                prop_mat_z%v_scr,                                         &
                                v_i,                                                      &
                                'n',it)
          Call Matrix_Vector_Multiply( prop_mat_d%triangle_hamiltonian,                   &
                                       prop_mat_z%v_scr,                                  &
                                       prop_mat_z%local_scratch)    
          Call Triangular_Solve(prop_mat_d%upper,                                         &
                                prop_mat_z%local_scratch,                                 &
                                v_o,                                                      &
                                't',it)
      END IF
  ELSE  ! This is for the standard case either with or without packed matrices
      IF (packed_matrices) THEN 
          Call Column_Packed_Symmetric_Matrix_on_Vector (                                 &
                             v_i,                                                         &
                             v_o,                                                         &
                             mat_var(3)%mat_d%packed_columns,                             &
                             mat_var(3)%non_zero_columns,                                 &
                             mat_var(3)%row_index,                                        &
                             mat_var(3)%mat_d%matrix_diagonal )
      END IF
      IF (matrix_type == 'triangle') THEN
          Call Matrix_Vector_Multiply( prop_mat_d%triangle_hamiltonian,                   &
                                       v_i,                                               &
                                       v_o)
      ELSE IF (matrix_type == 'fedvr') THEN
          Call finite_element_m_v(v_i,v_o)

      END IF
  END IF
END SUBROUTINE Generalized_Lanczos_Matmul_z
!***********************************************************************
!***********************************************************************
!**begin prologue     Generalized_Lanczos_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Perform the Lanczos recursion at a given timestep.
!**description        This routine does the actual work required at each timestep
!**                   
!**references
!**routines called    
!**end prologue       Generalized_Lanczos_d
!***********************************************************************
  SUBROUTINE Generalized_Lanczos_d(v_i,v_o)
  IMPLICIT NONE
  TYPE(REAL_PROP_VAR)                               :: prop_mat_d
  REAL(idp), DIMENSION(:)                           :: v_i
  REAL(idp), DIMENSION(:)                           :: v_o
  CHARACTER (LEN=5)                                 :: itoc
  INTEGER                                           :: it
  REAL(idp)                                         :: ddot
  CHARACTER(LEN=8)                                  :: kewrd
!
  null_vec = .false.
  convergence = 'unconverged'
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
  prop_mat_d%lanczos_tri=0
!
!
  local_time=secnds(0.0)
!
  v_i(:) = v_i(:) /sqrt(ddot(n3d,v_i(:),1,v_i(:),1))  ! First Lanczos vector
  prop_mat_d%vec(:,0) = v_i(:)                                 ! Copy into vector
  total_time(1) = secnds(0.0) - local_time + total_time(1)
!
!
  DO it=0, maxvec-1       ! Lanczos iterations
!
     WRITE(iout,3) it
!
!            Lets look at the input vector
!
     IF(log_iterative(10).or.log_iterative(2)) THEN
        call Print_Matrix(type_real_vector,prop_mat_d%vec(:,it),title='Lanczos Vector iteration = '//itoc(it))
     END IF
!
     local_time=secnds(0.0)
!
!    
!
     Call Generalized_Lanczos_Matmul(prop_mat_d%vec(:,it), prop_mat_d%h_vec(:), it) ! The crucial matrix vector multiply step
!
     total_time(2) = secnds(0.0) - local_time + total_time(2)
     IF(log_iterative(10).or.log_iterative(3)) THEN
        call Print_Matrix(type_real_vector,prop_mat_d%h_vec,title='H_on_Lanczos_Vector iteration = '//itoc(it))
     END IF
!
!             Compute the value of a(it) as the
!                               -T   -1
!                  a(i) = <V_i|U  H U  |V_i> 
!
     local_time=secnds(0.0)
     a(it) =  ddot(n3d , prop_mat_d%vec(:,it) , 1 , prop_mat_d%h_vec(:) , 1 )
!     IF(log_iterative(10).or.log_iterative(1)) THEN
        Write(iout,*) 'Lanczos Alpha iteration = '//itoc(it), a(it)
!     END IF
!
!
     Call Generalized_Lanczos_Hamiltonian(eig,eigen_vectors,sub_diagonal,                    &
                              eig_previous,eigen_vectors_previous,it) ! Compute eigenvalues and eigenvectors of the 
!                                                                     ! projected Hamiltonian
     total_time(3) = secnds(0.0) - local_time + total_time(3)
!
!
     local_time=secnds(0.0)
     call Generalized_Lanczos_Convergence_Test(v_o,eigen_vectors,prop_mat_d%lanczos_tri,it) ! Test convergence
     total_time(4) = secnds(0.0) - local_time + total_time(4)
!
!
     local_time = secnds(0.0)
     IF(convergence == 'unconverged') THEN ! Not converged generate new vector and continue
!
        Call Generalized_Lanczos_Vec(prop_mat_d%vec,prop_mat_d%h_vec,it)
     ELSE IF( convergence == 'maximal') THEN ! If the number of iterations has reached the size of the
                                             ! original matrix without breakdown, We have to be converged even if we
                                             ! have not reached the convergence criterion on the eignvalue.
                                             ! Perhaps we have set that criterion too tight for the numerics but in
                                             ! in any event we are finished and need to see what is going on.
!
        WRITE(iout,4) eig(0)
        norm = 1.d0 / sqrt(ddot(n3d,v_o(:),1,v_o(:),1))
        v_i(:) = norm * v_o(:)
        total_time(8) = secnds(0.0) - local_time + total_time(8)
        total_number_of_iterations = total_number_of_iterations + 1
        return
!
     ELSE IF(convergence == 'converged') THEN ! We actually have converged
!
        WRITE(iout,5) it  
        norm = 1.d0 / sqrt(ddot(n3d,v_o(:),1,v_o(:),1))
        v_i(:) = norm * v_o(:)
        total_number_of_iterations = total_number_of_iterations + 1
        total_time(8) = secnds(0.0) - local_time + total_time(8)
        return
    END IF
    total_number_of_iterations = total_number_of_iterations + 1
  END DO
  write(iout,6) maxvec
  norm = 1.d0 / sqrt(ddot(n3d,v_o(:),1,v_o(:),1))
  v_i(:) = norm * v_o(:)
  total_number_of_iterations = total_number_of_iterations + 1
  total_time(8) = secnds(0.0) - local_time + total_time(8)
!
1 FORMAT('***********************************************'  &
         '*************************') 
2 FORMAT(/,20X,'Beginning Lanczos Iterations At Zero')
3 FORMAT(/,20x,'Lanczos Iteration = ',i4)
4 FORMAT(/,20x,'Number of Iterations Maximal.  Energy = ', e15.8)
5 FORMAT(/,20x,'Convergence After ',i5, ' Iterations')
6 FORMAT(/,20x,'No Convergence After ',i5, 'Iterations')
END SUBROUTINE Generalized_Lanczos_d
!***********************************************************************
!***********************************************************************
!**begin prologue     Generalized_Lanczos_Hamiltonian
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Diagonalize tridiagonal matrix
!***                  
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Generalized_Lanczos_Hamiltonian
!***********************************************************************
  Subroutine Generalized_Lanczos_Hamiltonian(eig,eigen_vectors,sub_diagonal,eig_previous,eigen_vectors_previous,it)
  IMPLICIT NONE
  REAL(idp), DIMENSION(0:maxvec)                    :: eig
  REAL(idp), DIMENSION(0:maxvec)                    :: eig_previous
  REAL(idp), DIMENSION(0:maxvec)                    :: sub_diagonal  
  REAL(idp), DIMENSION(0:maxvec,0:maxvec)           :: eigen_vectors  
  REAL(idp), DIMENSION(0:maxvec,0:maxvec)           :: eigen_vectors_previous  
  INTEGER                                           :: it
  CHARACTER(LEN=4)                                  :: itoc
!
!               
!
!            Lets get the eigenvalues of the tridiagonal matix.
!
  IF(log_iterative(10).or.log_iterative(5)) THEN
     call Print_Matrix(type_real_vector,a(0:it),title='Diagonal Elements Lanczos H iteration = '//itoc(it) )
  END IF
  IF ( it > 0) THEN
!       write(iout,*) eig(0:it-1)
       sub_diagonal(0:it-1) = b(1:it)
       IF(log_iterative(10).or.log_iterative(5)) THEN
           call Print_Matrix(type_real_vector,sub_diagonal(0:it-1),        &
                             title = 'Off-Diagonal Element Lanczos H iteration = '//itoc(it))
       END IF
       eig_previous(0:it-1) = eig(0:it-1)
!       write(iout,*) eig_previous(0:it-1)
       eigen_vectors_previous(0:it-1,0:it-1) = eigen_vectors(0:it-1,0:it-1)
  END IF
  eig(0:it) = a(0:it)
  call dstev('v',it+1,eig,sub_diagonal,eigen_vectors,maxvec+1,rwork,info)
  IF ( it == 0 ) THEN
     eig_previous(0:0) = eig(0:0)
!       write(iout,*) eig_previous(0)
       eigen_vectors_previous(0:0,0:0) = eigen_vectors(0:0,0:0)
  END IF
  eig_old = eig(0)
  call Print_Matrix(type_real_vector,eig(0:it),title='Lanczos eigenvalues iteration = '//itoc(it))
  IF(log_iterative(10).or.log_iterative(7)) THEN
     call Print_Matrix(type_real_matrix,eigen_vectors(0:it,0:it), it+1,it+1,                       &
                       title='Eigenvectors in lanczos basis iteration = '//itoc(it))
  END IF
!
END SUBROUTINE Generalized_Lanczos_Hamiltonian
!***********************************************************************
!***********************************************************************
  SUBROUTINE Generalized_Lanczos_Convergence_Test_d(v_o,eigen_vectors,lanczos_tri,it)
  IMPLICIT NONE
  TYPE(REAL_PROP_VAR)                               :: prop_mat_d
  REAL(idp), DIMENSION(:)                           :: v_o
  REAL(idp), DIMENSION(0:maxvec,0:maxvec)           :: eigen_vectors
  REAL(idp), DIMENSION(0:maxvec)                    :: lanczos_tri
  INTEGER                                           :: i, it
  INTEGER                                           :: t_interval
  REAL(idp)                                         :: ddot
  CHARACTER(LEN=4)                                  :: itoc
!
!                     This is the exponentiation step
!  
!             Compute the projection of the time-dependent state on the initial vector,
!             The initial vector is the first Lanczos vector.  So, this just scales
!             the first component of each eigenvector by an exponential.
!
!             C(t) = Sum_{lambda} Exp(-lambda * deltat ) d_{lambda}(0) C_{lambda}
!                                                     +
!             Compute the d_{lambda}(0) = [C_{lambda}]  S C(0)
!
  v_o(1:it+1) = eigen_vectors(0,0:it)
  Call Exponential_on_Vector(prop_mat_d%work,v_o,eig,eigen_vectors,it)    
  IF ( it + 1 == n3d) THEN  
       lanczos_tri(0:it) = lanczos_tri(0:it) - prop_mat_d%work(0:it)
       wfn_tst = sqrt ( ddot( it+1 , lanczos_tri , 1 ,lanczos_tri , 1 ) )
       write(iout,1) wfn_tst
       call ebcx(v_o,n3d,prop_mat_d%vec,n3d,eigen_vectors,maxvec+1,n3d,it+1,1)
       convergence = 'maximal'
       return
  END IF
  t_end = t_start + deltat
  IF( it+1 == maxvec) THEN
      lanczos_convergence = 'unconverged'
      t_interval=0
      save_deltat = deltat
      DO i = 1, maximum_number_of_time_subintervals
         Write(iout,2)
         lanczos_tri(0:it) = 0.d0
         deltat=.5d0*deltat
         v_o(1:it) = eigen_vectors_previous(0,0:it-1)
         Call Exponential_on_Vector(prop_mat_d%work,v_o,eig_previous,eigen_vectors_previous,it-1)    
         lanczos_tri(0:it-1) = prop_mat_d%work(0:it-1)                  
!         write(iout,*) work_d(0:it-1)
         v_o(1:it+1) = eigen_vectors(0,0:it)
         Call Exponential_on_Vector(prop_mat_d%work,v_o,eig,eigen_vectors,it)    
         lanczos_tri(0:it) = lanczos_tri(0:it) - prop_mat_d%work(0:it)                  
!         write(iout,*) work_d(0:it)
         wfn_tst = sqrt ( ddot( it+1 , lanczos_tri , 1 ,lanczos_tri , 1 ) )
         t_interval = t_interval + 1 
         t_end = t_start + deltat
         write(iout,3) t_interval, wfn_tst, deltat, t_end
         IF (wfn_tst <= cnverg ) THEN
             lanczos_convergence = 'converged'
             convergence = 'converged'
             call ebcx(v_o,n3d,prop_mat_d%vec,n3d,eigen_vectors,maxvec+1,n3d,it+1,1)
             exit
         END IF
      END DO
      write(iout,*) 'Final Value of Lowest Eigenvalue = ',eig(0)
      eig_old = eig(0)
      deltat = save_deltat
      write(iout,4) t_interval, wfn_tst, t_end
      call ebcx(v_o,n3d,prop_mat_d%vec,n3d,eigen_vectors,maxvec+1,n3d,it+1,1)
      return
  END IF
!
!             work_d contains the exponentiated vector in the Lanczos basis.  local scratch has the effect of the exponential
!             on the starting vector.

!             At this point, if we were converged, we would transform to the
!             original basis and return.  However, if we are not converged
!             we do not have to do that.  So, lets test convergence.
!
!
!  IF( it+1 == maxvec) THEN
!
!             We have maximized the number of iterations, we are totally finished.  
!             Take the lowest eigenvalue and eigenvector, transform to the original basis and return,
!               
!       convergence = 'maximal'
!       title='Final Value of Lowest Eigenvalue'
!       call prntfmn(title,eig,1,1,maxvec+1,1,iout,'e')
!       eig_old = eig(0)
!       IF(log_iterative(10).or.log_iterative(9)) THEN     
!          call ebcx(v_o,n3d,vec_d,n3d,work_d,maxvec+1,n3d,it+1,1)
!          title='Final Value of Exp(-H*t) on Initial Vector'
!          call prntfmn(title,v_o,n3d,1,n3d,1,iout,'e')
!       END IF 
!       call ebcx(v_o,n3d,vec_d,n3d,eigen_vectors,maxvec+1,n3d,it+1,1)
!       lanczos_tri(0:it) = lanczos_tri(0:it) - work_d(0:it)
!       wfn_tst = sqrt ( ddot( it+1 , lanczos_tri , 1 ,lanczos_tri , 1 ) )
!       write(iout,1) wfn_tst
!       lanczos_tri(0:it) = work_d(0:it)
!       RETURN
!  END IF
  IF (it == 0) THEN
       lanczos_tri(it) = prop_mat_d%work(it)
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          call Print_Matrix(type_real_vector,lanczos_tri(0:it),title='lanczos_tri iteration = '//itoc(it))
       END IF
  ELSE
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          call Print_Matrix(type_real_vector,lanczos_tri(0:it),title='lanczos_tri iteration = '//itoc(it))
          call Print_Matrix(type_real_vector,prop_mat_d%work(0:it),title='work iteration = '//itoc(it))
       END IF
       lanczos_tri(0:it) = lanczos_tri(0:it) - prop_mat_d%work(0:it)
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          call Print_Matrix(type_real_vector,lanczos_tri(0:it),title='lanczos_tri iteration = '//itoc(it))
       END IF
  END IF
  wfn_tst = sqrt ( ddot( it+1 , lanczos_tri , 1 ,lanczos_tri , 1 ) )
  write(iout,1) wfn_tst
  lanczos_tri(0:it) = prop_mat_d%work(0:it)
!
!             We are below the maximal number of iterations.  Test convergence.
!
  IF(wfn_tst <= cnverg) THEN
!
!                  
!    If converged return eigenvalue and eigenvector
!    in original basis.
!
       convergence = 'converged'
       write(iout,*) 'Converged Value of Lowest Eigenvalue = ', eig(0)
!
!    Transform to original basis
!    C_alpha(t) = Sum vec_d(alpha,i) * A(i)
!    C(t) = Sum C_alpha(t) * q_alpha
!
     eig_old = eig(0)
     call ebcx(v_o,n3d,prop_mat_d%vec,n3d,eigen_vectors,maxvec+1,n3d,it+1,1)
     IF(log_iterative(10).or.log_iterative(9)) THEN     
         call Print_Matrix(type_real_vector,v_o,title='Final Solution')
     END IF
     RETURN
  END IF
1 FORMAT(/,10x,'RMS = sqrt(< ( v_{i} - v_{i-1}) |( v_{i} - v_{i-1}) >) = ',e15.8)
2 FORMAT(/,25x,'Reducing the Time Step')
3 FORMAT(10x,'Sub-Iteration = ',i3,'  RMS = sqrt(< ( v_{i} - v_{i-1}) |( v_{i} - v_{i-1}) >) = ',e15.8,  &
               /,10x,'New Step Size = ',e15.8,' Final Time = ',e15.8)
4 FORMAT(10x,'Iteration = ',i3,'  RMS = sqrt(< ( v_{i} - v_{i-1}) |( v_{i} - v_{i-1}) >) = ',e15.8,      &
             /,10x,'Converged',' Final Time = ',e15.8)
END SUBROUTINE Generalized_Lanczos_Convergence_Test_d
!***********************************************************************
!***********************************************************************
!**begin prologue     Generalized_Lanczos_Vec_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Compute next Lanczos vector.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Lanczos_Vec_d
!***********************************************************************
  SUBROUTINE Generalized_Lanczos_Vec_d(lanczos_basis,h_vec,it)
  IMPLICIT NONE
  REAL(idp), DIMENSION(n3d,0:maxvec)                :: lanczos_basis
  REAL(idp), DIMENSION(n3d)                         :: h_vec
  INTEGER                                           :: i
  INTEGER                                           :: it
  INTEGER                                           :: low
  REAL(idp)                                         :: ddot
  CHARACTER(LEN=4)                                  :: itoc
!
!            Note that it runs from 0 to maxvec
!            The starting vector is assumed to be S normalized.
!
!
!
!          Compute next vector at it + 1
!
!                 First compute
!
!    b_{i+1} |V_{i+1}> = [ H - a_{i}  ] |V_{i}> -   b_{i}|V_{i-1}>
!    and store it in V_{i+1}
!
     lanczos_basis(:,it+1) = h_vec(:)  - a(it) * lanczos_basis(:,it)
!
     IF( it > 0 ) THEN
!
         lanczos_basis(:,it+1) = lanczos_basis(:,it+1) - b(it) * lanczos_basis(:,it-1)
!
         IF(log_iterative(10)) THEN
            call Print_Matrix(type_real_vector,lanczos_basis(:,it+1),title='before reorthogonalization')
         END IF
!
!        Perform a reorthogonalization to prevent linear dependence.
!
         local_time=secnds(0.0)
         low=0
         IF(orthogonalize =='double_schmidt') THEN
             low = it - 1
         ELSE IF(orthogonalize =='full') THEN
             low = 0
         ELSE IF(orthogonalize =='partial') THEN
            low = min(it_min,it-1)
         END IF
         DO i=low, it
            lanczos_basis(:,it+1) = lanczos_basis(:,it+1)                        &
                                 -                                               &
                       ddot(n3d,lanczos_basis(:,i),1,lanczos_basis(:,it+1),1)    &
                                 * lanczos_basis(:,i)
         END DO
!
         total_time(5) = secnds(0.0) - local_time + total_time(5)
         IF(log_iterative(10)) THEN
            call Print_Matrix(type_real_vector,lanczos_basis(:,it+1),title='after reorthogonalization')
         END IF
     END IF
!
!        To compute the next beta, we use the formula,
!
!                      beta(it+1)  = Sqrt( <lanczos_basis(:,it+1) | lanczos_basis(:,it+1) > )
!
!        This comes from the normalization requirement on the new vector.
!
     local_time=secnds(0.0)
!
!
     b(it+1) = SQRT( ddot(n3d , lanczos_basis(:,it+1) , 1 , lanczos_basis(:,it+1) , 1) ) 
!
!     IF(log_iterative(10).or.log_iterative(2)) THEN
        Write(iout,*) 'Lanczos Beta iteration = '//itoc(it), b(it+1)
!     END IF
!
     local_time=secnds(0.0)
!
!          Form the new basis function
!
     IF( b(it+1) >  eps ) THEN
         lanczos_basis(:,it+1) =  lanczos_basis(:,it+1) / b(it+1)
         total_time(7) = secnds(0.0) - local_time + total_time(7)
     ELSE IF( b(it+1) <= eps ) THEN
!
!
!        Test for breakdown.  It occurs if b(i) becomes zero.
!
         local_time = secnds(0.0)
         null_vec=.true.
!
!        Compute a random vector and then gram-schmidt it to the previous lanczos
!        vectors.
!
         write(iout,1)
         call random_number(lanczos_basis(:,it+1))
         DO i=0,it
            lanczos_basis(:,it+1) = lanczos_basis(:,it+1)                              &
                                         -                                             &
                             ddot( n3d,lanczos_basis(:,i),1,lanczos_basis(:,it+1),1 )  &
                                 * lanczos_basis(:,i)
         END DO
         b(it+1) = SQRT( ddot(n3d , lanczos_basis(:,it+1) , 1 , lanczos_basis(:,it+1) , 1) ) 
         lanczos_basis(:,it+1) = lanczos_basis(:,it+1) / b(it+1)
         IF(log_iterative(10).or.log_iterative(2)) THEN
            Write(iout,*) 'New Lanczos Beta iteration = '//itoc(it+1), b(it+1)
         END IF
         total_time(6) = secnds(0.0) - local_time + total_time(6)
         b(it+1) = 0.d0
     END IF
!
!
1 FORMAT(/,10x,'Lanczos Breakdown.  Generate New Vector and Continue')
END SUBROUTINE Generalized_Lanczos_Vec_d
!***********************************************************************
!***********************************************************************
!**begin prologue     Generalized_Lanczos_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Perform the Lanczos recursion at a given timestep.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Generalized_Lanczos_z
!***********************************************************************
  SUBROUTINE Generalized_Lanczos_z(v_i,v_o)
  IMPLICIT NONE
  TYPE(COMPLEX_PROP_VAR)                :: prop_mat_z
  COMPLEX(idp), DIMENSION(:)            :: v_i
  COMPLEX(idp), DIMENSION(:)            :: v_o
  CHARACTER (LEN=5)                     :: itoc
  INTEGER                               :: it
  COMPLEX(idp)                          :: cdotc
  CHARACTER(LEN=8)                      :: kewrd
!
  eps=1.d-04
  null_vec = .false.
  convergence = 'unconverged'
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
  prop_mat_z%lanczos_tri=0
!
!
  local_time=secnds(0.0)
  v_i(:) = v_i(:) / sqrt( cdotc(n3d, v_i, 1, v_i, 1))
  prop_mat_z%vec(:,0) = v_i(:)
  total_time(1) = secnds(0.0) - local_time + total_time(1)
!
!           Begin Lanczos iterations
!
  DO it=0, maxvec-1
!
     WRITE(iout,3) it
!
!
!            Lets look at the input vector
!
     IF(log_iterative(10).or.log_iterative(2)) THEN
        call Print_Matrix(type_complex_vector,prop_mat_z%vec(:,it),title='Lanczos Vector iteration = '//itoc(it))
     END IF
!
     local_time=secnds(0.0)
!
     Call Generalized_Lanczos_Matmul(prop_mat_z%vec(:,it), prop_mat_z%h_vec(:), it)
!
     total_time(2) = secnds(0.0) - local_time + total_time(2)
     IF(log_iterative(10).or.log_iterative(3)) THEN
        call Print_Matrix(type_complex_vector,prop_mat_z%h_vec,title='H_on_Lanczos_Vector iteration = '//itoc(it))
     END IF
!
!             Compute the value of a(it) as the
!                  a(i) = <V_i|h|V_i> 
!
     local_time=secnds(0.0)
     a(it) =  cdotc(n3d , prop_mat_z%vec(:,it) , 1 , prop_mat_z%h_vec(:) , 1 )
     IF(log_iterative(10).or.log_iterative(1)) THEN
        Write(iout,*) 'Lanczos Alpha iteration = '//itoc(it), a(it)
     END IF
!
!             Compute eigenvalues and eigenvectors of the projected
!             Hamiltonian
!
     Call Generalized_Lanczos_Hamiltonian(eig,eigen_vectors,sub_diagonal,                    &
                              eig_previous,eigen_vectors_previous,it)
     total_time(3) = secnds(0.0) - local_time + total_time(3)
!
!
!            Test Convergence
!
     local_time=secnds(0.0)
     call Generalized_Lanczos_Convergence_Test(v_o,eigen_vectors,prop_mat_z%lanczos_tri,it)
     total_time(4) = secnds(0.0) - local_time + total_time(4)
!
    local_time = secnds(0.0)
    IF( convergence == 'unconverged') THEN
!
!           We are not converged and need to generate a new vector and continue.
!
        Call Generalized_Lanczos_Vec(prop_mat_z%vec,prop_mat_z%h_vec,it)
!
    ELSE IF( convergence == 'maximal') THEN
!
!           If the number of iterations has reached the size of the
!           original matrix without breakdown, we have to be converged.
!
!
             WRITE(iout,4) eig(0)
             norm = 1.d0/sqrt(cdotc(n3d,v_o(:),1,v_o(:),1))
             v_i(:) = norm * v_o(:)
             IF(compute_energy) THEN 
                Call Generalized_Lanczos_Matmul(v_i(:), prop_mat_z%h_vec(:), it)
                eig(1) = cdotc(n3d,v_i(:),1,prop_mat_z%h_vec(:),1)
                write(iout,*) 'Expectation Value of Hamiltonian = ', eig(1)
             END IF
             total_time(8) = secnds(0.0) - local_time + total_time(8)
             total_number_of_iterations = total_number_of_iterations + it + 1
             return
!
!            We are actually converged to the desired tolerance.
!
    ELSE IF(convergence == 'converged') THEN
            WRITE(iout,5) it  
            norm = 1.d0/sqrt(cdotc(n3d,v_o(:),1,v_o(:),1))
            v_i(:) = norm * v_o(:)
            total_number_of_iterations = total_number_of_iterations + it + 1
            return 
            IF(compute_energy) THEN 
               Call Generalized_Lanczos_Matmul(v_i(:), prop_mat_z%h_vec(:), it)
            END IF
            eig(1) = cdotc(n3d,v_i(:),1,prop_mat_z%h_vec(:),1)
            write(iout,*) 'Expectation Value of Hamiltonian = ', eig(1)
    END IF
    total_time(8) = secnds(0.0) - local_time + total_time(8)
    total_number_of_iterations = total_number_of_iterations + it + 1
  END DO
  write(iout,6) maxvec
  norm = 1.d0/sqrt(cdotc(n3d,v_o(:),1,v_o(:),1))
  v_i(:) = norm * v_o(:)
  total_time(8) = secnds(0.0) - local_time + total_time(8)
1 FORMAT('***********************************************'  &
         '*************************') 
2 FORMAT(/,5X,'Beginning Lanczos Iterations Beginning At Zero')
3 FORMAT(/,10x,'Lanczos Iteration = ',i4)
4 FORMAT(/,20x,'Number of Iterations Maximal')
5 FORMAT(/,20x,'Convergence After ',i5, ' Iterations')
6 FORMAT(/,20x,'No Convergence After ',i5, 'Iterations')
END SUBROUTINE Generalized_Lanczos_z
!***********************************************************************
!***********************************************************************
!**begin prologue     Generalized_Lanczos_Convergence_Test_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Perform the Lanczos recursion at a given timestep.
!***                  b_{i+1} |V_{i+1}> = [ H - a_{i} ] |V_{i}> - b_{i} |V_{i-1}>
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Generalized_Lanczos_Convergence_Test_z
!***********************************************************************
  SUBROUTINE Generalized_Lanczos_Convergence_Test_z(v_o,eigen_vectors,lanczos_tri,it)
  IMPLICIT NONE
  TYPE(COMPLEX_PROP_VAR)                               :: prop_mat_z
  REAL(idp), DIMENSION(0:maxvec,0:maxvec)              :: eigen_vectors
  COMPLEX(idp), DIMENSION(0:maxvec)                    :: lanczos_tri
  COMPLEX(idp), DIMENSION(:)                           :: v_o
  COMPLEX(idp)                                         :: temp
  INTEGER                                              :: i, it
  INTEGER                                              :: t_interval
  COMPLEX(idp)                                         :: cdotc
  CHARACTER(LEN=4)                                     :: itoc
!
!               
  v_o(1:it+1) = eigen_vectors(0,0:it)
  Call Exponential_on_Vector(prop_mat_z%work,v_o,eig,eigen_vectors,it)
  IF ( it + 1 == n3d) THEN  
       lanczos_tri(0:it) = lanczos_tri(0:it) - prop_mat_z%work(0:it)
       wfn_tst = sqrt ( cdotc( it+1 , lanczos_tri , 1 ,lanczos_tri , 1 ) )
       call cebcx(v_o,n3d,prop_mat_z%vec,n3d,eigen_vectors,maxvec+1,n3d,it+1,1)
       convergence = 'maximal'
       return
  END IF
!
!          At this point, if we were converged, we would transform to the
!          original basis and return.  However, if we are not converged
!          we do not have to do that.  So, lets test convergence.
  t_end = t_start + deltat
  IF( it+1 == maxvec) THEN
      lanczos_convergence = 'unconverged'
      t_interval=0
      save_deltat = deltat
      DO i = 1, maximum_number_of_time_subintervals
         Write(iout,2)
         lanczos_tri(0:it) = 0.d0
         deltat=.5d0*deltat
         v_o(1:it) = eigen_vectors_previous(0,0:it-1)
         Call Exponential_on_Vector(prop_mat_z%work,v_o,eig_previous,eigen_vectors_previous,it-1)    
         lanczos_tri(0:it-1) = prop_mat_z%work(0:it-1)                  
!         write(iout,*) prop_mat_z%work(0:it-1)
         v_o(1:it+1) = eigen_vectors(0,0:it)
         Call Exponential_on_Vector(prop_mat_z%work,v_o,eig,eigen_vectors,it)    
         lanczos_tri(0:it) = lanczos_tri(0:it) - prop_mat_z%work(0:it)                  
!         write(iout,*) prop_mat_z%work(0:it)
         wfn_tst = sqrt ( cdotc( it+1 , lanczos_tri , 1 ,lanczos_tri , 1 ) )
         t_interval = t_interval + 1 
         t_end = t_start + deltat
         write(iout,3) t_interval, wfn_tst, deltat, t_end
         IF (wfn_tst <= cnverg ) THEN
             lanczos_convergence = 'converged'
             convergence = 'converged'
             call cebcx(v_o,n3d,prop_mat_z%vec,n3d,prop_mat_z%work,maxvec+1,n3d,it+1,1)
             exit
         END IF
      END DO
      deltat = save_deltat
      write(iout,4) t_interval, wfn_tst, t_end
      call cebcx(v_o,n3d,prop_mat_z%vec,n3d,prop_mat_z%work,maxvec+1,n3d,it+1,1)
      call Print_Matrix(type_complex_vector,v_o,title='Final Value of Exp(-i*H*t) on Initial Vector Final')
      IF(log_iterative(10).or.log_iterative(9)) THEN     
         call Print_Matrix(type_complex_vector,v_o,title='Final Solution')
      END IF 
      RETURN
  END IF
  IF (it == 0) THEN
       lanczos_tri(it) = prop_mat_z%work(it)
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          call Print_Matrix(type_complex_vector,lanczos_tri(0:it),title='lanczos_tri iteration = '//itoc(it))
       END IF
  ELSE
      IF(log_iterative(10).or.log_iterative(8)) THEN     
         call Print_Matrix(type_complex_vector,lanczos_tri(0:it),title='lanczos_tri iteration = '//itoc(it))
         call Print_Matrix(type_complex_vector,prop_mat_z%work(0:it),title='work iteration = '//itoc(it))
      END IF
      lanczos_tri(0:it) = lanczos_tri(0:it) - prop_mat_z%work(0:it)
      IF(log_iterative(10).or.log_iterative(8)) THEN     
         call Print_Matrix(type_complex_vector,lanczos_tri(0:it),title='lanczos_tri iteration = '//itoc(it))
      END IF
  END IF
  wfn_tst = sqrt ( cdotc( it+1 , lanczos_tri , 1 ,lanczos_tri , 1 ) )
  write(iout,1) wfn_tst
  lanczos_tri(0:it) = prop_mat_z%work(0:it)
!
  IF(wfn_tst <= cnverg) THEN
!
!                  
!    If converged return eigenvalue and eigenvector
!    in original basis.
!
       convergence = 'converged'
!
!    Transform to original basis
!    C_alpha(t) = Sum vec_d(alpha,i) * A(i)
!    C(t) = Sum C_alpha(t) * q_alpha
!
     call cebcx(v_o,n3d,prop_mat_z%vec,n3d,prop_mat_z%work,maxvec+1,n3d,it+1,1)
     IF(log_iterative(10).or.log_iterative(9)) THEN     
        call Print_Matrix(type_complex_vector,v_o,title='Final Solution')
     END IF
     RETURN
  END IF
1 FORMAT(/,10x,'RMS = sqrt(< ( v_{i} - v_{i-1}) |( v_{i} - v_{i-1}) >) = ',e15.8)
2 FORMAT(/,25x,'Reducing the Time Step')
3 FORMAT(10x,'Sub-Iteration = ',i3,'  RMS = sqrt(< ( v_{i} - v_{i-1}) |( v_{i} - v_{i-1}) >) = ',e15.8,  &
               /,10x,'New Step Size = ',e15.8,' Final Time = ',e15.8)
4 FORMAT(10x,'Iteration = ',i3,'  RMS = sqrt(< ( v_{i} - v_{i-1}) |( v_{i} - v_{i-1}) >) = ',e15.8,      &
             /,10x,'Converged',' Final Time = ',e15.8)
END SUBROUTINE Generalized_Lanczos_Convergence_Test_z
!***********************************************************************
!***********************************************************************
!**begin prologue     Generalized_Lanczos_Vec_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Perform the Lanczos recursion at a given timestep.
!***                  b_{i+1} |V_{i+1}> = [ H - a_{i} ] |V_{i}> - b_{i} |V_{i-1}>
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Generalized_Lanczos_Vec_z
!***********************************************************************
  SUBROUTINE Generalized_Lanczos_Vec_z(lanczos_basis,h_vec,it)
  IMPLICIT NONE
  COMPLEX(idp), DIMENSION(n3d,0:maxvec)            :: lanczos_basis
  COMPLEX(idp), DIMENSION(n3d)                     :: h_vec
  INTEGER                                          :: i
  INTEGER                                          :: it
  INTEGER                                          :: low
  COMPLEX(idp)                                     :: cdotc
  CHARACTER(LEN=4)                                 :: itoc
  REAL(idp)                                        :: t_ran
!
!            Note that it runs from 0 to maxvec
!            The starting vector is assumed to be normalized.
!
!
!
!          Compute next vector at it + 1
!
!                 First compute
!
!    b_{i+1} |V_{i+1}> = [ H - a_{i}  ] |V_{i}> -   b_{i}|V_{i-1}>
!    and store it in V_{i+1}
!
     lanczos_basis(:,it+1) = h_vec(:)  - a(it) * lanczos_basis(:,it)
!
     IF( it > 0 ) THEN
!
!
         lanczos_basis(:,it+1) = lanczos_basis(:,it+1) - b(it) * lanczos_basis(:,it-1)
!
         IF(log_iterative(10)) THEN

            call Print_Matrix(type_complex_vector,lanczos_basis(:,it+1),title='before reorthogonalization')
         END IF
!
!        Perform a reorthogonalization to prevent linear dependence.
!
         local_time=secnds(0.0)
         low=0
         IF(orthogonalize=='double_schmidt') THEN
             low = it - 1
         ELSE IF(orthogonalize=='full') THEN
             low = 0
         ELSE IF(orthogonalize=='partial') THEN
            low = min(it_min,it-1)
         END IF
         DO i=low, it
            lanczos_basis(:,it+1) = lanczos_basis(:,it+1)                           &
                                 -                                                  &
                       cdotc(n3d,lanczos_basis(:,i),1,lanczos_basis(:,it+1),1 )     &
                                 * lanczos_basis(:,i)
         END DO
!
         total_time(5) = secnds(0.0) - local_time + total_time(5)
         IF(log_iterative(10)) THEN
            call Print_Matrix(type_complex_vector,lanczos_basis(:,it+1),title='after reorthogonalization')
         END IF
     END IF
!
!        To compute the next beta, we use the formula,
!
!                     2          
!                 beta  = <v_o|v_o) 
!
!        This comes from the S normalization requirement on the new vector.
!
     local_time=secnds(0.0)
!
!         Compute b(it+1)  b(i) = Sqrt ( <V_i|V_i> )
!
     b(it+1) = SQRT( cdotc(n3d , lanczos_basis(:,it+1) , 1 , lanczos_basis(:,it+1) , 1) ) 
!
!     IF(log_iterative(10).or.log_iterative(2)) THEN
        write(iout,*) 'Lanczos Beta iteration = '//itoc(it), b(it+1)
!     END IF
!
!          Form the new basis function
!
     IF( b(it+1) >  eps ) THEN
         lanczos_basis(:,it+1) =  lanczos_basis(:,it+1) / b(it+1)
         total_time(7) = secnds(0.0) - local_time + total_time(7)
     ELSE IF( b(it+1) <= eps ) THEN
!
!        Test for breakdown.  It occurs if b(i) becomes zero.
!
         local_time = secnds(0.0)
         null_vec=.true.
         write(iout,1)
         DO i=1,n3d
            call random_number(t_ran)
            lanczos_basis(i,it+1) = t_ran
         END DO
         DO i=0,it
            lanczos_basis(:,it+1) = lanczos_basis(:,it+1)                               &
                                          -                                             &
                             cdotc( n3d,lanczos_basis(:,i),1,lanczos_basis(:,it+1),1 )  &
                                 * lanczos_basis(:,i)
         END DO
         b(it+1) = SQRT( cdotc(n3d , lanczos_basis(:,it+1) , 1 , lanczos_basis(:,it+1) , 1) ) 
         lanczos_basis(:,it+1) = lanczos_basis(:,it+1)  / b(it+1)
         IF(log_iterative(10).or.log_iterative(2)) THEN
            Write(iout,*) 'New Lanczos Beta iteration = '//itoc(it+1), b(it+1)
         END IF
         total_time(6) = secnds(0.0) - local_time + total_time(6)
         b(it+1) = 0.d0
     END IF
!
!
1 FORMAT(/,10x,'Lanczos Breakdown.  Generate New Vector and Continue')
END SUBROUTINE Generalized_Lanczos_Vec_z
!***********************************************************************
!***********************************************************************
           END MODULE Generalized_Lanczos_Module
!***********************************************************************
!***********************************************************************
