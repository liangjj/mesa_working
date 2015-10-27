!***********************************************************************
! Propagation_Module
!**begin prologue     Propagation_Module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Arnoldi, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains all of the major subroutines to propagate
!***                  a wavefunction in time using the Lanczos algorithm.
!***                  Explicit interfaces are used to allow
!***                  a transparent use of generic subroutines which work
!***                  for both real and complex vectors.  This feature
!***                  permits a single code to be used for both real and
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
!--------------------------------------------------------------------------------------------
!***                  Note, there are three cases treated in this package.  The standard
!***                  symmetric problem, where the overlap matrix is unity, the generalized
!***                  symmetric problem but explicitly transformed to a standard problem
!***                  using the Cholesky decomposition and what is termed the preconditioned
!***                  version where the Lanczos iterations are performed using the Cholesky
!***                  decomposition at each iteration by solving a set of triangular linear
!***                  equations.
!--------------------------------------------------------------------------------------------
!**references
!**modules needed     See USE statements below
!**comments           In this portable version I have disabled all unnecessary
!**                   writing to files.  The original Fortran is commented out.
!**                   In addition, there is no option to compute the autocorrelation
!**                   function as this would require reading and manipulating the
!**                   initial state wavefunction from a file.
!**end prologue       Propagation_Module
!***********************************************************************
!***********************************************************************
                           MODULE Propagation_Module
                           USE initial_state_module
                           USE Atomic_Matrices,                                          &    
                                                              ONLY:                      &
                                                              atomic_eigen_values,       &
                                                              state_matrix_size,         &
                                                              final_state_matrix_size,   &
                                                              final_total_matrix_size,   &
                                                              state_mat,                 &
                                                              dipole_mat
!
                           USE Global_Time_Propagation_Module,                           &
                                                              ONLY:                      &
                                                              dipole_matrices,           &
                                                              non_orth,                  &
                                                              time_dependent_potential,  &
                                                              packed_matrices,           &
                                                              use_atomic_symmetry,       &
                                                              scratch_matrix,            &
                                                              state_file_name,           &
                                                              packed_file_name,          &
                                                              file_key

                           USE Generalized_Lanczos_Module
                           USE Iterative_Global,                                         &
                                                              ONLY:                      &
                                                              local_scratch_d,           &             
                                                              local_scratch_z   
                           USE Preconditioned_Lanczos_Module
                           USE auto_correlation_module
                           USE moment_module
                           USE Data_Subroutines_Module,       ONLY : len
                           USE input_output
                           USE FEDVR_Shared, ONLY : file_loc
!***********************************************************************
!
                          INTERFACE Propagator
                   MODULE PROCEDURE Propagator_d,                          &
                                    Propagator_z
                          END INTERFACE Propagator
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Iterative_Data
!**begin prologue     Iterative_Data
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, iterative, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Read the data required for the iteration process.
!***                  Parameters specifying the maximum number of iteration vectors, 
!***                  convergence criterion and overlap tolerance for vectors kept.
!***description       reads in parameters for iteration.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Iterative_Data
  SUBROUTINE Iterative_Data
  IMPLICIT NONE
  LOGICAL                       :: dollar, logkey
  CHARACTER (LEN=80)            :: chrkey
  REAL*8                        :: fpkey
  INTEGER                       :: intkey 
!
! Read in the data
!
  IF ( dollar('$iterative_data',card,cpass,inp) ) then
      iterative_method=chrkey(card,'iterative_method','standard_lanczos',' ')
      print_iterative='print=iterative='//print_iterative
      print_iterative(10)=chrkey(card,'print=iterative=',print_iterative(10),' ')
      IF(print_iterative(10) == 'all') THEN
         CALL setprn(log_iterative,10)
      ELSE
         CALL setlog(log_iterative,print_iterative,card,9)
      END IF
!
!     When cnverg is reached the iterations stop.
!
      cnverg=fpkey(card,'convergence',1.d-08,' ')
!     
!     If the overlap matrix of an vector gets too small
!     it means linear dependence may be creeping in.  The vector
!     needs to be discarded and this may signal a breakdown in the
!     process.  The user needs to be careful about the quality of the results.
!
      thresh=fpkey(card,'overlap_tolerance',1.d-08,' ')
      orthogonalize = chrkey(card,'orthogonalization_procedure','none',' ')
      it_min = intkey(card,'number_of_iterates',0,' ')
      non_orth=logkey(card,'non_orthogonal_basis',.false.,' ')
      type=chrkey(card,'type_matrix','general',' ')
      compute_energy = logkey(card,'compute_energy',.false.,' ')
!
!     Maximum number of iterations and vectors to keep.
!
      maxit=intkey(card,'maximum_number_of_iterations',n3d,' ')
      maxvec=intkey(card,'maximum_number_of_vectors',n3d,' ')
      maxit=MIN(maxit,n3d)
      maxvec=MIN(maxvec,n3d)
!
!     Number of trial vectors and a now obsolete restart capability.
!     The restart never worked.
!
      WRITE(iout,1) thresh, cnverg, maxit, maxvec, type, non_orth, orthogonalize
  ELSE
      write(iout,2)
      stop
  END IF
1 FORMAT(/,15X,'iterative diagonalization information',/,&
         /,5X,'overlap tolerance                  = ',e15.8, &
         /,5X,'convergence criterion              = ',e15.8, &
         /,5X,'maximum number of iterations       = ',i6,    &
         /,5X,'maximum number of vectors          = ',i6,    &
         /,5X,'type of Hamiltonian                = ',a16,   &
         /,5X,'non_orthogonal basis               = ',l1,    &
         /,5X,'orthogonalization procedure        = ',a16)
2 FORMAT(/,5x,'error in diagonalization card section')
END SUBROUTINE Iterative_Data
!**********************************************************************
!**********************************************************************
!                            Driving Routine
!**********************************************************************
!***********************************************************************
! Propagation_Driver
!**begin prologue     Propagation_Driver
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains all of the major subroutines to propagate
!***                  a wavefunction in time using the Arnoldi/Lanczos
!***                  algorithm.  Explicit interfaces are used to allow
!***                  a transparent use of generic subroutines which work
!***                  for both real and complex vectors.  This feature
!***                  permits a single code to be used for both real and
!***                  imaginary time propagation.
!***description       First, we get an initial vector.  There are a number of choices
!***                  possible. Given the starting vector, we then call a 
!***                  routine to allocate storage.  The allocation depends on
!***                  the method chosen, that is whether it is a standard Lanczos,
!***                  generalized or preconditioned Lanczos.  After that the
!***                  actual routine that does the work is called.  In that routine
!***                  a number of ierations are performed until the time propagated solution
!***                  satisfies a fixed accuracy criterion.  The stopping
!***                  criterion used is based on how accurate the vector
!***                  has converged at that timestep.  In imaginary time, 
!***                  the number of timesteps is controlled by the convergence
!***                  of the eigenvalue of interest.  The converged vector 
!***                  provides the starting vector for the
!***                  next time step.
!**references
!**modules needed     See USE statements below
!**comments           In this portable version I have disabled all unnecessary
!**                   writing to files.  The original Fortran is commented out.
!**                   In addition, there is no option to compute the autocorrelation
!**                   function as this would require reading and manipulating the
!**                   initial state wavefunction from a file.
!**end prologue       Propagation_Driver
!***********************************************************************
!***********************************************************************
  SUBROUTINE Propagation_Driver
  IMPLICIT NONE
  INTEGER                       :: i 
  lwork = 5*(maxvec+1)
!
! Allocate memory for the main arrays.              
!                                                                         
! Get Data on Initial State.  Many possibilities exist.  
!
  CALL Initial_State_Data    ! Default is unit vector.
!
  IF ( type_calculation == 'real_time' ) THEN
       ALLOCATE( psi_z(1:n3d), v_scr_z(1:n3d) )
       IF ( i0stat == 'unperturbed_state_vector'                        &
                           .or.                                         &
            i0stat == 'perturbed_state_vector') THEN                   
            ALLOCATE(itemp(n3d,spdim),etemp(n3d))        
            CALL Initial_Vector(psi_z)          ! Fill up the initial vector
            DEALLOCATE(itemp,etemp)        
       ELSE IF( i0stat == 'superpose') THEN                           
            ALLOCATE(c_loc(spdim))
            DO i=1,spdim
               ALLOCATE(c_loc(i)%c(nphy(i)),c_loc(i)%phi(nphy(i)),      &
                        c_loc(i)%l(nphy(i)))
            END DO
            CALL Initial_Vector(psi_z)          ! Fill up the initial vector
            DO i=1,spdim
               DEALLOCATE(c_loc(i)%c,c_loc(i)%phi,c_loc(i)%l)
            END DO
            DEALLOCATE(c_loc)
       ELSE
            CALL Initial_Vector(psi_z)          ! Fill up the initial vector
       END IF
       Call Propagator_Setup                    ! Set up arrays
       Call Propagator ( psi_z, v_scr_z )       ! Do the heavy lifting
       DEALLOCATE( psi_z, v_scr_z )
  ELSE IF ( type_calculation == 'imaginary_time' ) THEN
       ALLOCATE( psi_d(1:n3d), v_scr_d(1:n3d) )
       IF ( i0stat == 'unperturbed_state_vector'                       &
                          .or.                                         &
            i0stat == 'perturbed_state_vector') THEN                   
            ALLOCATE(itemp(n3d,spdim),etemp(n3d))        
            CALL Initial_Vector(psi_d)          ! Fill up the initial vector
            DEALLOCATE(itemp,etemp)        
       ELSE IF( i0stat == 'superpose') THEN                           
            ALLOCATE(c_loc(spdim))
            DO i=1,spdim
               ALLOCATE(c_loc(i)%c(nphy(i)),c_loc(i)%phi(nphy(i)),     &
                        c_loc(i)%l(nphy(i)))
            END DO
            CALL Initial_Vector(psi_d)          ! Fill up the initial vector
            DO i=1,spdim
               DEALLOCATE(c_loc(i)%c,c_loc(i)%phi,c_loc(i)%l)
            END DO
            DEALLOCATE(c_loc)
       ELSE
            CALL Initial_Vector(psi_d)          ! Fill up the initial vector
       END IF
       Call Propagator_Setup                    ! Set up arrays
       Call Propagator (psi_d, v_scr_d )        ! Do the hevy lifting
       DEALLOCATE( psi_d, v_scr_d )
  ELSE
       call lnkerr('Quit.  Bad Time Keyword')
  END IF
 END SUBROUTINE Propagation_Driver
!**********************************************************************
!**********************************************************************
!deck propagator_setup
!**begin prologue     propagator_setup
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, propagate
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Time dependent Schroedinger equation
!**                   using iterative method with finite difference or
!**                   dvr space representation.  This subroutine is for
!**                   real vectors.
!***description       This routine allocates the storage that is needed in the
!***                  propagation with the exception of the two arrays that hold
!***                  the wavefunction and a needed scratch vector.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       propagator_setup
  SUBROUTINE propagator_setup
  IMPLICIT NONE
  INTEGER                             :: i
  INTEGER                             :: L
  INTEGER                             :: count
  CHARACTER (LEN=8)                   :: itoc
  WRITE(iout,1)
  IF ( type_calculation == 'imaginary_time' ) THEN
!                  
!                        All these arrays are needed
!
       ALLOCATE( vec_d(n3d,0:maxvec), h_vec_d(n3d) )
       ALLOCATE( eig(0:maxvec), sub_diagonal(0:maxvec),                                    &
                 eigen_vectors(0:maxvec,0:maxvec), lanczos_tri_d(0:maxvec),                &
                 work_d(0:lwork), rwork(0:lwork), eig_previous(0:maxvec),                  &
                 eigen_vectors_previous(0:maxvec,0:maxvec), a(0:maxvec), b(0:maxvec))
       IF(iterative_method == 'standard_lanczos') THEN
          Write(iout,2)
       ELSE IF(iterative_method == 'generalized_lanczos') THEN  ! An additional array is required
                                                                ! for the generalized Lanczos
          WRITE(iout,3)
          ALLOCATE(local_scratch_d(n3d))
       ELSE IF(iterative_method == 'preconditioned_lanczos') THEN  ! An additional matrix is required
                                                                   ! for the preconditioned Lanczos
          WRITE(iout,4)
          ALLOCATE( s_vec_d(n3d,0:maxvec) )
       END IF
!
  ELSE IF(type_calculation == 'real_time') THEN
       ALLOCATE( vec_z(n3d,0:maxvec), h_vec_z(n3d))
       ALLOCATE( eig(0:maxvec), sub_diagonal(0:maxvec),                                    &
                 eigen_vectors(0:maxvec,0:maxvec), lanczos_tri_z(0:maxvec),                &
                 rwork(0:lwork), work_z(0:lwork), eig_previous(0:maxvec),                  &
                 eigen_vectors_previous(0:maxvec,0:maxvec), a(0:maxvec), b(0:maxvec) )
       IF(iterative_method == 'standard_lanczos') THEN
          Write(iout,2)
          IF (time_dependent_potential == .true. .and. dipole_matrices == .true. ) THEN
              IF (use_atomic_symmetry == 'on' ) THEN   ! default is off and it is unclear
                                                       ! that the on mode actually works yet
!
!                 All of the stuff below is in place when one has matrices on a state-by-state-basis
!                 These must already exist.  There is an option in the code to take the original
!                 matrices and format them this way.  If that option is invoked, the code will reformat
!                 and then restart with the formatted matrices and go on.
!
                  write(iout,*) '*** Reopening File to Read State Matrices ***'
                  write(iout,*) '    *** Reading Needed Variables From State Matrices ***'
                  file_loc=File_Directory(4)(1:len_dir(4))//'/'//'state_matrix_file'
                  Call pakstr(file_loc,len(1))
                  Call IOsys('open state_matrices as old',0,0,0,file_loc(1:len(1)))
                  Call iosys('read integer L_Max from state_matrices',1,L_Max,0,' ')
                  ALLOCATE( state_matrix_size(0:L_Max),                                        &
                            final_state_matrix_size(0:L_Max)                         ,         &
                            dipole_mat(0:L_Max) )
                  Call IOsys('read integer state_matrix_size from state_matrices',             &
                              L_max+1,state_matrix_size, 0 , ' ' )       
                  Call IOsys('close state_matrices',0,0,0,' ')          
!
!                 Now we take the packed matrix file, reopen it and fill it with the required
!                 dipole arrays
!
                  write(iout,*) '*** Reopening File to Read Packed Matrices ***'
                  write(iout,*) '    *** Reading Needed Variables From Packed Matrices ***'
                  file_loc=File_Directory(4)(1:len_dir(4))//'/'//'packed_matrix_file'
                  Call pakstr(file_loc,len(1))
                  Call IOsys('open packed_matrices as old',0,0,0,file_loc(1:len(1)))
                  Call IOsys('read integer final_state_matrix_size from packed_matrices',      &
                              L_max+1,final_state_matrix_size, 0 , ' ' )       
                  Call IOsys('read integer final_total_matrix_size from packed_matrices',      &
                              1,final_total_matrix_size, 0 ,' ' )       
                  Call IOsys('read real "atomic_eigenvalues" from packed_matrices',            &
                              final_total_matrix_size,atomic_eigen_values,0,' ')
                  DO L = 0, L_max
                     file_key = 'L_'//itoc(L)
                     Call IOsys('read integer "number_non_zero_c_dipole_matrix_'//file_key//   &
                                '_elements" from packed_matrices',1,                           &
                                dipole_mat(L)%dipole_number,0,' ')
                     ALLOCATE (                                                                &
                               dipole_mat(L)%dipole_non_zero_columns                           &
                                                          ( final_state_matrix_size(L+1) ),    &
                               dipole_mat(L)%dipole_row_index                                  &
                                                          ( dipole_mat(L)%dipole_number ),     &
                               dipole_mat(L)%dipole_packed_columns                             &
                                                          ( dipole_mat(L)%dipole_number) )                    
                     Call Read_and_Write_Column_Packed_Matrices (                              &
                                                        dipole_mat(L)%dipole_packed_columns,   &
                                                        dipole_mat(L)%dipole_non_zero_columns, &
                                                        dipole_mat(L)%dipole_row_index,        &
                                                        'read',                                &
                                                        'c_dipole_matrix_'//file_key,          &
                                                        dipole_mat(L)%dipole_number )
                  END DO
              END IF
          END IF   
       ELSE IF(iterative_method == 'generalized_lanczos') THEN   
              ALLOCATE(local_scratch_z(n3d))     ! allocate some additional storage
              WRITE(iout,3)
       ELSE IF(iterative_method == 'preconditioned_lanczos') THEN
              WRITE(iout,4)
              ALLOCATE( s_vec_z(n3d,0:maxvec) )  ! allocate some additional storage
       END IF
  END IF
!
1 FORMAT('***********************************************'                                     &
         '*************************')
2 FORMAT(/,20x,'Begin Standard Lanczos Propagation',/10x,                                      &
               'Orthonormal Basis, No Preconditioning')
3 FORMAT(/,20x,'Begin Generalized Lanczos Propagation',/10x,                                   &
               'Generalized EIgenvalue Problem Transformed to Standard Form')
4 FORMAT(/,20x,'Begin Preconditioned Lanczos Propagation',/10x,                                &
               'Generalized Eigenvalue Problem Direct Solution' )
END SUBROUTINE propagator_setup
!**********************************************************************
!**********************************************************************
!deck Propagator_d
!**begin prologue     Propagator_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, propagate
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Time dependent Schroedinger equation
!**                   using iterative method with finite difference or
!**                   dvr space representation.  This subroutine is for
!**                   real vectors.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Propagator_d
  SUBROUTINE Propagator_d(input_vector, scratch_vector)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                :: input_vector                         
  REAL*8, DIMENSION(:)                :: scratch_vector                       
  INTEGER                             :: t
  INTEGER                             :: i
  CHARACTER (LEN=8)                   :: kewrd
  WRITE(iout,1)
  ALLOCATE(total_time(10))
  total_time = 0.0
  t_start = t_init
  t_end = t_start
  eig_tst = 0.d0
  total_number_of_iterations = 0
  DO  t=1,ntreg   ! Begin the propagation 
!
!     Propagate to the next time using the short time formula,
!\begin{equation}
!    \Xi(x,y,z,t_{i-1} + {\delta}t) = \big [ exp( - i H(x,y,z,t_{i-1}) {\delta}t )
!                                -1 \big ] \Psi_{0}(x,y,z,t_{i-1})
!                   \nonumber
!\end{equation}
!
!     This is done in the Lanczos code by performing an effective reduction
!     of the Hamiltonian to diagonal form using the value of the wavefunction
!     from the previous step as an initial guess.
!
       IF( iterative_method == 'generalized_lanczos') THEN
           CALL generalized_lanczos(input_vector,scratch_vector)
       ELSE IF( iterative_method == 'preconditioned_lanczos') THEN
           CALL preconditioned_lanczos(input_vector,scratch_vector)
       END IF
       IF (log_main(8)) THEN
           call plot_psi(input_vector,t,t_end)
       END IF
       write(iout,2) t, norm, eig_old
       eig_tst = abs(eig_tst-eig_old)
       IF(convergence == 'maximal') THEN  ! We cannot continue as we have exceeded the maximum number of iterations
                                          ! We cannot get a more accurate eigenvalue.
!      
          WRITE(iout,1)
          WRITE(iout,3)
          WRITE(iout,1)
          EXIT                           ! Got to quit
       END IF
       IF( eig_tst <= cnverg ) THEN      ! Bingo, we converged
!
!     End the propagation
!
           WRITE(iout,1)
           WRITE(iout,4)
           WRITE(iout,1)
           EXIT
       ELSE
           eig_tst = eig_old            ! Replace the test value of the eigenvalue by the current value
                                        ! and coninue
       END IF
       t_start=t_end
  END DO
!
  WRITE(iout,*) 'Writing Final Wavefunction to Disk'
  kewrd='real'
  REWIND(20)
  WRITE(20) kewrd
  WRITE(20) input_vector
  CLOSE(20)
!    End the propagation
!
  WRITE(iout,1)
  WRITE(iout,5) total_time(1:8), total_number_of_iterations
  WRITE(iout,1)
!
!    Deallocate Memory in preparation for normalexit
!
  DEALLOCATE( vec_d, h_vec_d )
  DEALLOCATE( eig, sub_diagonal, eigen_vectors, lanczos_tri_d, &
              rwork, work_d, eig_previous,eigen_vectors_previous, a, b )
  IF(iterative_method == 'standard_lanczos') THEN
!
     IF (packed_matrices) THEN
         DEALLOCATE( mat_var(3)%non_zero_columns,      &                       
                     mat_var(3)%row_index,             &
                     mat_var(3)%packed_columns_d,      &
                     mat_var(3)%matrix_diagonal_d )
         DEALLOCATE(mat_var)
     END IF
!
  ELSE IF(iterative_method == 'preconditioned_lanczos') THEN
     IF (packed_matrices ) THEN
         DEALLOCATE(mat_var(1)%non_zero_columns,       &
                    mat_var(1)%row_index,              &
                    mat_var(1)%packed_columns_d,       &
                    mat_var(1)%matrix_diagonal_d )
         DEALLOCATE(mat_var(2)%non_zero_columns,       &
                    mat_var(2)%row_index,              &
                    mat_var(2)%packed_columns_d )
         DEALLOCATE(mat_var(3)%non_zero_columns,       &
                    mat_var(3)%row_index,              &
                    mat_var(3)%packed_columns_d,       &
                    mat_var(3)%matrix_diagonal_d )
         DEALLOCATE(mat_var)
     ELSE
         DEALLOCATE(upper_d)
     END IF
!
  END IF
!
  DEALLOCATE(total_time)
  RETURN
!
1 FORMAT('***********************************************'             &
         '*************************')
2 FORMAT(/,5x,'Time Step = ',i6,/,10x,'Normalization = ',e15.8,1x,     &
              'Energy    = ',e15.8)
3 FORMAT(/,20x,'Iterations Maximal : Stop')
4 FORMAT(/,20x,'Converged : End Propagation')
5 FORMAT(/,20x,'Detailed Timing',/,20x,                                 &
               'Initialization             = ',e15.8,/,20x,             &
               'Matrix Vector Multiply     = ',e15.8,/,20x,             &
               'Hamiltonian Construction   = ',e15.8,/,20x,             &
               'Convergence Test           = ',e15.8,/,20x,             &
               'Schmidt Process            = ',e15.8,/,20x,             &
               'Lanczos Breakdown          = ',e15.8,/,20x,             &
               'Linear System Solves       = ',e15.8,/,20x,             &
               'Cleanup                    = ',e15.8,/,20x,             &
               'Total Number of Iterations = ',i10)              
END SUBROUTINE Propagator_d
!***********************************************************************
!***********************************************************************
!deck Propagator_z
!**begin prologue     Propagator_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagate
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Time dependent Schroedinger equation
!**                   using iterative method with finite difference or
!**                   dvr space representation.  This routine is for complex
!**                   vectors.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       Propagator_z
  SUBROUTINE Propagator_z(input_vector, scratch_vector)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)            :: input_vector                         
  COMPLEX*16, DIMENSION(:)            :: scratch_vector                       
  INTEGER                             :: t
  INTEGER                             :: i
  CHARACTER (LEN=8)                   :: kewrd
!
!                   Will not repeat the comments as they are basically the same as in the real routine
!
  ALLOCATE(total_time(10))
  total_time = 0.0
  total_number_of_iterations = 0
  t_start = t_init
  t_end = t_start
  DO  t=1,ntreg
!
!        Propagate to the next time using the short time formula,
!\begin{equation}
!    \Xi(x,y,z,t_{i-1} + {\delta}t) = \big [ exp( - i H(x,y,z,t_{i-1}) {\delta}t )
!                                -1 \big ] \Psi_{0}(x,y,z,t_{i-1})
!                   \nonumber
!\end{equation}
!        This is done in the Lanczos code by performing an effective reduction
!        of the Hamiltonian to diagonal form using the value of the wavefunction
!        from the previous step as an initial guess.
!
      IF( iterative_method == 'generalized_lanczos') THEN
          IF (use_atomic_symmetry == 'on') THEN
!               Call atomic_lanczos(input_vector, scratch_vector)
          ELSE
               CALL generalized_lanczos(input_vector,scratch_vector)
          END IF
      ELSE IF( iterative_method == 'preconditioned_lanczos') THEN
         CALL preconditioned_lanczos(input_vector,scratch_vector)
      END IF
      input_vector = input_vector/sqrt(norm)
      IF (log_main(8)) THEN
          call plot_psi(input_vector,t,t_end)
      END IF
      IF(i0stat == 'gaussian_pulse') THEN
         call calculate_moment(input_vector,t1)
      END IF
      write(iout,2) t, norm
!  
      t_start=t_end
  END DO
  WRITE(iout,*) 'Writing Final Wavefunction to Disk'
  kewrd='complex'
  REWIND(20)
  WRITE(20) kewrd
  WRITE(20) input_vector
  CLOSE(20)
  WRITE(iout,1)
  WRITE(iout,3)
  WRITE(iout,1)
!
!    End the propagation
!
!    Deallocate Memory
!
  DEALLOCATE( vec_z, h_vec_z )
  DEALLOCATE( eig, sub_diagonal, eigen_vectors, lanczos_tri_z,   &
              work_z, eig_previous,eigen_vectors_previous, a, b )
  IF(iterative_method == 'generalized_lanczos') THEN
!
     IF (packed_matrices) THEN
         DEALLOCATE( mat_var(3)%non_zero_columns,      &                       
                     mat_var(3)%row_index,             &
                     mat_var(3)%packed_columns_z,      &
                     mat_var(3)%matrix_diagonal_z )
         DEALLOCATE(mat_var)
     END IF
!
  ELSE IF(iterative_method == 'preconditioned_lanczos') THEN
     IF (packed_matrices ) THEN
         DEALLOCATE(mat_var(1)%non_zero_columns,       &
                    mat_var(1)%row_index,              &
                    mat_var(1)%packed_columns_z,       &
                    mat_var(1)%matrix_diagonal_z )
         DEALLOCATE(mat_var(2)%non_zero_columns,       &
                    mat_var(2)%row_index,              &
                    mat_var(2)%packed_columns_z )
         DEALLOCATE(mat_var(3)%non_zero_columns,       &
                    mat_var(3)%row_index,              &
                    mat_var(3)%packed_columns_z,       &
                    mat_var(3)%matrix_diagonal_z )
         DEALLOCATE(mat_var)
     ELSE
         DEALLOCATE(upper_z)
     END IF
  END IF
  WRITE(iout,4) total_time(1:8), total_number_of_iterations
  DEALLOCATE(total_time)
  RETURN
1 FORMAT('***********************************************'             &
         '*************************')
2 FORMAT(/,5x,'Time Step = ',i6,/,10x,'Normalization = ',e15.8)
3 FORMAT(/,20x,'Converged: End Propagation')
4 FORMAT(/,20x,'Detailed Timing',/,20x,                                 &
               'Initialization             = ',e15.8,/,20x,             &
               'Matrix Vector Multiply     = ',e15.8,/,20x,             &
               'Hamiltonian Construction   = ',e15.8,/,20x,             &
               'Convergence Test           = ',e15.8,/,20x,             &
               'Schmidt Process            = ',e15.8,/,20x,             &
               'Linear System Solves       = ',e15.8,/,20x,             &
               'Lanczos Breakdown          = ',e15.8,/,20x,             &
               'Cleanup                    = ',e15.8,/,20x,             &
               'Total Number of Iterations = ',i10)              
END SUBROUTINE Propagator_z
!***********************************************************************
!***********************************************************************
END MODULE Propagation_Module
!***********************************************************************
!***********************************************************************
