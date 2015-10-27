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
!***                  a wavefunction in time using the Iterative or Arnoldi
!***                  algorithms.  Explicit interfaces are used to allow
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
                           USE Exact_Propagator_Module
                           USE Lanczos_Module
                           USE Arnoldi_Module
                           USE auto_correlation_module
                           USE moment_module
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
                           INTERFACE propagation
             MODULE PROCEDURE propagation_d,                       &
                              propagation_z
                       END INTERFACE propagation
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
!**purpose            set up the parameters specifying the maximum
!**                   number of iteration vectors, convergence criterion and
!***                  overlap tolerance for vectors kept.
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
      iterative_method=chrkey(card,'iterative_method','lanczos',' ')
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
!***description       Given a starting vector, a number of ierations
!***                  are performed until the time propagated solution
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
! Get Data on Initial State
!
  CALL Initial_State_Data
!
  IF ( type_calculation == 'real_time' ) THEN
       ALLOCATE( psi_z(1:n3d),v_scr_z(1:n3d) ,v_tot(n3d))
       IF ( i0stat == 'unperturbed_state_vector'                       &
                           .or.                                        &
            i0stat == 'perturbed_state_vector') THEN                   
            ALLOCATE(itemp(n3d,spdim),etemp(n3d))        
            CALL Initial_Vector(psi_z)
            DEALLOCATE(itemp,etemp)        
       ELSE IF( i0stat == 'superpose') THEN                           
            ALLOCATE(c_loc(spdim))
            DO i=1,spdim
               ALLOCATE(c_loc(i)%c(nphy(i)),c_loc(i)%phi(nphy(i)),     &
                        c_loc(i)%l(nphy(i)))
            END DO
            CALL Initial_Vector(psi_z)
            DO i=1,spdim
               DEALLOCATE(c_loc(i)%c,c_loc(i)%phi,c_loc(i)%l)
            END DO
            DEALLOCATE(c_loc)
       ELSE
            CALL Initial_Vector(psi_z)
       END IF
       CALL Propagation ( psi_z, v_scr_z )
       DEALLOCATE( psi_z, v_scr_z, v_tot )
  ELSE IF(type_calculation == 'imaginary_time' ) THEN
       ALLOCATE( psi_d(1:n3d),v_scr_d(1:n3d) ,v_tot(n3d))
       IF ( i0stat == 'unperturbed_state_vector'                       &
                          .or.                                         &
            i0stat == 'perturbed_state_vector') THEN                   
            ALLOCATE(itemp(n3d,spdim),etemp(n3d))        
            CALL Initial_Vector(psi_d)
            DEALLOCATE(itemp,etemp)        
       ELSE IF( i0stat == 'superpose') THEN                           
            ALLOCATE(c_loc(spdim))
            DO i=1,spdim
               ALLOCATE(c_loc(i)%c(nphy(i)),c_loc(i)%phi(nphy(i)),     &
                        c_loc(i)%l(nphy(i)))
            END DO
            CALL Initial_Vector(psi_d)
            DO i=1,spdim
               DEALLOCATE(c_loc(i)%c,c_loc(i)%phi,c_loc(i)%l)
            END DO
            DEALLOCATE(c_loc)
       ELSE
            CALL Initial_Vector(psi_d)
       END IF
       CALL Propagation( psi_d, v_scr_d )
       DEALLOCATE( psi_d, v_scr_d, v_tot )
  ELSE
       call lnkerr('Quit.  Bad Time Keyword')
  END IF
1 FORMAT(/5x,'Allocating storage for psi_z,v_scr_z,v_tot')
END SUBROUTINE Propagation_Driver
!**********************************************************************
!**********************************************************************
!deck propagation_d
!**begin prologue     propagation_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Time dependent Schroedinger equation
!**                   using iterative method with finite difference or
!**                   dvr space representation.  This subroutine is for
!**                   real vectors.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       propagation_d
  SUBROUTINE propagation_d(wave_function,scratch_vector)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                :: wave_function
  REAL*8, DIMENSION(:)                :: scratch_vector
  CHARACTER (LEN=2)                   :: itoc
  INTEGER                             :: i, t
  INTEGER                             :: lenth, len
  REAL*8                              :: t_cur
  REAL*8                              :: ddot
  WRITE(iout,1)
  ALLOCATE(total_time(10))
  total_time = 0.0
  ALLOCATE( vec_d(n3d,0:maxvec), h_vec_d(n3d), eig(0:maxvec),          &  
            sub_diagonal(0:maxvec), eigen_vectors(0:maxvec,0:maxvec),  &
            lanczos_tri_d(0:maxvec), work_d(0:lwork), rwork(0:lwork))
  IF(iterative_method == 'lanczos') THEN
     WRITE(iout,2)
!
!    Allocate memory for the lanczos subroutines.              
!
     ALLOCATE( a(0:maxvec), b(0:maxvec))
  ELSE IF(iterative_method == 'preconditioned_lanczos') THEN
     WRITE(iout,3)
!
!    Allocate memory for the preconditioned lanczos subroutines.              
!
     ALLOCATE( a(0:maxvec), b(0:maxvec), s_vec_d(n3d,0:maxvec) )
  ELSE IF(iterative_method == 'arnoldi') THEN
     WRITE(iout,4)
!
!    Allocate memory for the Arnoldi subroutines.              
!
     ALLOCATE( h_mat_work_d(0:maxvec,0:maxvec) )
  ELSE IF(iterative_method == 'preconditioned_arnoldi') THEN
     WRITE(iout,5)
     call lnkerr('quit. Not yet implimented')
  ELSE
     call lnkerr('Quit.  Bad Method Keyword')
  END IF
  WRITE(iout,1)
  t_cur=t_init
  eig_tst =0.d0
  REWIND 99  
  READ(99) wave_function
  total_number_of_iterations = 0
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=lenth(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=lenth(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=t_cur
      t1=t0+deltat
!
!        Calculate the time dependent perturbation.
!        It consists of a space and a time part.
      v_tot = 0.d0
      IF (.not.no_pot) THEN
          CALL v_tim
          CALL pert
      END IF
!  
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
     IF( iterative_method == 'lanczos') THEN
         CALL lanczos(wave_function,scratch_vector)
     ELSE IF( iterative_method == 'preconditioned_lanczos') THEN
         CALL preconditioned_lanczos(wave_function,scratch_vector)
     ELSE IF (iterative_method == 'arnoldi') THEN
         CALL arnoldi(wave_function,scratch_vector)
     ELSE IF( iterative_method == 'preconditioned_arnoldi') THEN
         Call lnkerr('quit. Not yet implimented') 
     END IF
     wave_function = wave_function/sqrt(norm)
     IF (log_main(8)) THEN
         call plot_psi(wave_function,t,t1)
     END IF
     write(iout,6) t, norm, eig_old
     eig_tst = abs(eig_tst-eig_old)
     IF(convergence == 'maximal') THEN

!       We do not need to continue as the maximal number of iterations have been achieved
!       and we cannot get a more accurate eigenvalue.
!       
        WRITE(iout,1)
        WRITE(iout,7)
        WRITE(iout,1)
        WRITE(iout,9) total_time(1:8), total_number_of_iterations
        WRITE(iout,1)
        RETURN
     END IF
     IF( eig_tst <= cnverg ) THEN
!
!        End the propagation
!
          WRITE(iout,1)
          WRITE(iout,8)
          WRITE(iout,1)
          WRITE(iout,9) total_time(1:8), total_number_of_iterations
          WRITE(iout,1)
          RETURN
     ELSE
            eig_tst = eig_old
!            WRITE(iout,10)
!            WRITE(iout,9) total_time(1:8), total_number_of_iterations
     END IF
!     END IF
  t_cur=t1
  END DO
!
!    End the propagation
!
!    Deallocate Memory
!
!
  DEALLOCATE( vec_d, h_vec_d, eig, sub_diagonal, eigen_vectors, lanczos_tri_d, rwork, work_d )

  IF(iterative_method == 'lanczos') THEN
!
     DEALLOCATE( a, b )
     IF (packed) THEN
         DEALLOCATE( mat_var(3)%non_zero_columns,      &                       
                     mat_var(3)%row_index,             &
                     mat_var(3)%packed_columns_d,      &
                     mat_var(3)%matrix_diagonal_d )
         DEALLOCATE(mat_var)
     END IF
!
  ELSE IF(iterative_method == 'preconditioned_lanczos') THEN
     DEALLOCATE( a, b, s_vec_d)
     IF (packed ) THEN
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
         DEALLOCATE(mat_var(4)%non_zero_columns,       &
                    mat_var(4)%row_index,              &
                    mat_var(4)%packed_columns_d,       &
                    mat_var(4)%matrix_diagonal_d )
         DEALLOCATE(mat_var)
     ELSE
         DEALLOCATE(upper_d, lower_d )
     END IF
!
  ELSE IF(iterative_method == 'arnoldi') THEN
!
     DEALLOCATE( h_mat_work_d )
  ELSE IF(iterative_method == 'preconditioned_arnoldi') THEN
  END IF
!
  DEALLOCATE(total_time)
!
1 FORMAT('***********************************************'             &
         '*************************')
2 FORMAT(/,20x,'Begin Lanczos Propagation')
3 FORMAT(/,20x,'Begin Preconditioned Lanczos Propagation')
4 FORMAT(/,20x,'Begin Arnoldi Propagation')
5 FORMAT(/,20x,'Begin Preconditioned Arnoldi Propagation')
6 FORMAT(/,5x,'Time Step = ',i6,/,10x,'Normalization = ',e15.8,1x,     &
              'Energy    = ',e15.8)
7 FORMAT(/,20x,'Iterations Maximal : End Propagation')
8 FORMAT(/,20x,'Converged : End Propagation')
9 FORMAT(/,20x,'Detailed Timing',/,20x,                                 &
               'Initialization             = ',e15.8,/,20x,             &
               'Matrix Vector Multiply     = ',e15.8,/,20x,             &
               'Hamiltonian Construction   = ',e15.8,/,20x,             &
               'Convergence Test           = ',e15.8,/,20x,             &
               'Schmidt Process            = ',e15.8,/,20x,             &
               'Lanczos Breakdown          = ',e15.8,/,20x,             &
               'Linear System Solves       = ',e15.8,/,20x,             &
               'Cleanup                    = ',e15.8,/,20x,             &
               'total Number of Iterations = ',i10)              
10 FORMAT(/,20x,'Time Steps Exhausted : End Propagation')
END SUBROUTINE propagation_d
!***********************************************************************
!***********************************************************************
!deck propagation_z
!**begin prologue     propagation_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Time dependent Schroedinger equation
!**                   using iterative method with finite difference or
!**                   dvr space representation.  This routine is for complex
!**                   vectors.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       propagation_z
  SUBROUTINE propagation_z(wave_function,scratch_vector)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)            :: wave_function
  COMPLEX*16, DIMENSION(:)            :: scratch_vector
  CHARACTER (LEN=2)                   :: itoc
  INTEGER                             :: i, t
  INTEGER                             :: lenth, len
  REAL*8                              :: t_cur
  COMPLEX*16                          :: cdotc
  ALLOCATE(total_time(10))
  WRITE(iout,1)
  ALLOCATE( vec_z(n3d,0:maxvec), h_vec_z(n3d), eig(0:maxvec),          &
            sub_diagonal(0:maxvec), eigen_vectors(0:maxvec,0:maxvec),  &
            lanczos_tri_z(0:maxvec), rwork(0:lwork), work_z(0:lwork))
  IF(iterative_method == 'lanczos') THEN
     WRITE(iout,2)
!
!    Allocate memory for the lanczos subroutines.              
!
     ALLOCATE( a(0:maxvec), b(0:maxvec))
  ELSE IF(iterative_method == 'preconditioned_lanczos') THEN
     WRITE(iout,3)
!
!    Allocate memory for the preconditioned lanczos subroutines.              
!
     ALLOCATE( a(0:maxvec), b(0:maxvec), s_vec_z(n3d,0:maxvec) )
  ELSE IF(iterative_method == 'arnoldi') THEN
     WRITE(iout,4)
!
!    Allocate memory for the iterative subroutines.              
!
     ALLOCATE( h_mat_work_z(0:maxvec,0:maxvec) )
  ELSE IF(iterative_method == 'preconditioned_arnoldi') THEN
     WRITE(iout,5)
  ELSE
     call lnkerr('Quit.  Bad Method Keyword')
  END IF
  WRITE(iout,1)
  REWIND 99
  READ(99) wave_function
  t_cur=t_init
  DO  t=1,ntreg
      keywrd=itoc(t)
      LEN=lenth(keywrd)
      keywrd='t'//keywrd(1:LEN)
      LEN=lenth(keywrd)
      keywrd='$v0('//keywrd(1:LEN)//')'
      t0=t_cur
      t1=t0+deltat
  
!        Calculate the time dependent perturbation.
!        It consists of a space and a time part.

      v_tot = 0.d0
      IF (.not.no_pot) THEN
          CALL v_tim
          CALL pert
      END IF
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
      IF( iterative_method == 'lanczos') THEN
          CALL lanczos(wave_function,scratch_vector)
      ELSE IF( iterative_method == 'preconditioned_lanczos') THEN
         CALL preconditioned_lanczos(wave_function,scratch_vector)
      ELSE IF (iterative_method == 'arnoldi') THEN
          CALL arnoldi(wave_function,scratch_vector)
      ELSE IF( iterative_method == 'preconditioned_arnoldi') THEN
      END IF
      wave_function = wave_function/sqrt(norm)
      IF (log_main(8)) THEN
          call plot_psi(wave_function,t,t1)
      END IF
      IF(i0stat == 'gaussian_pulse') THEN
         call calculate_moment(wave_function,t1)
      END IF
      write(iout,6) t, norm
!  
      t_cur=t1
  END DO
  WRITE(iout,1)
  WRITE(iout,7)
  WRITE(iout,1)
!
!    End the propagation
!
!    Deallocate Memory
!
!
  DEALLOCATE( vec_z, h_vec_z, eig, sub_diagonal, eigen_vectors, lanczos_tri_z, work_z )
  IF(iterative_method == 'lanczos') THEN
!
     DEALLOCATE( a, b )
     IF (packed) THEN
         DEALLOCATE( mat_var(3)%non_zero_columns,      &                       
                     mat_var(3)%row_index,             &
                     mat_var(3)%packed_columns_z,      &
                     mat_var(3)%matrix_diagonal_z )
         DEALLOCATE(mat_var)
     END IF
!
  ELSE IF(iterative_method == 'preconditioned_lanczos') THEN
     DEALLOCATE( a, b, s_vec_z)
     IF (packed ) THEN
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
         DEALLOCATE(mat_var(4)%non_zero_columns,       &
                    mat_var(4)%row_index,              &
                    mat_var(4)%packed_columns_z,       &
                    mat_var(4)%matrix_diagonal_z )
         DEALLOCATE(mat_var)
     ELSE
         DEALLOCATE(upper_z, lower_z )
     END IF
  ELSE IF(iterative_method == 'arnoldi') THEN
!
     DEALLOCATE( h_mat_work_z )
  ELSE IF(iterative_method == 'preconditioned_arnoldi') THEN
  END IF
  WRITE(iout,8) total_time(1:8), total_number_of_iterations
  DEALLOCATE(total_time)
1 FORMAT('***********************************************'             &
         '*************************')
2 FORMAT(/,20x,'Begin Lanczos Propagation')
3 FORMAT(/,20x,'Begin Preconditioned Lanczos Propagation')
4 FORMAT(/,20x,'Begin Arnoldi Propagation')
5 FORMAT(/,20x,'Begin Preconditioned Arnoldi Propagation')
6 FORMAT(/,5x,'Time Step = ',i6,/,10x,'Normalization = ',e15.8)
7 FORMAT(/,20x,'Converged: End Propagation')
8 FORMAT(/,20x,'Detailed Timing',/,20x,                                 &
               'Initialization             = ',e15.8,/,20x,             &
               'Matrix Vector Multiply     = ',e15.8,/,20x,             &
               'Hamiltonian Construction   = ',e15.8,/,20x,             &
               'Convergence Test           = ',e15.8,/,20x,             &
               'Schmidt Process            = ',e15.8,/,20x,             &
               'Linear System Solves       = ',e15.8,/,20x,             &
               'Lanczos Breakdown          = ',e15.8,/,20x,             &
               'Cleanup                    = ',e15.8,/,20x,             &
               'total Number of Iterations = ',i10)              
END SUBROUTINE propagation_z
!***********************************************************************
!***********************************************************************
END MODULE Propagation_Module
!***********************************************************************
!***********************************************************************
