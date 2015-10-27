!***********************************************************************
! Lanczos_Module
!**begin prologue     Lanczos_Module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains all of the major subroutines to propagate
!***                  a wavefunction in time using the Lanczos
!***                  algorithm.  Explicit interfaces are used to allow
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
!**end prologue       Lanczos_Module
!***********************************************************************
!***********************************************************************
                           MODULE Lanczos_Module
                           USE full_matrix_vector_multiply_module
                           USE packed_matrix_vector_multiply_module
                           USE dvr_matrix_vector_multiply_module
                           USE Preconditioner_Module
                   INTEGER                            :: total_number_of_iterations
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                           INTERFACE lanczos
             MODULE PROCEDURE lanczos_d,                           &
                              lanczos_z
                       END INTERFACE lanczos
!
                           INTERFACE lanczos_vec
             MODULE PROCEDURE lanczos_vec_d,                       &
                              lanczos_vec_z
                       END INTERFACE lanczos_vec
!
                           INTERFACE lanczos_hamiltonian
             MODULE PROCEDURE lanczos_hamiltonian              
                       END INTERFACE lanczos_hamiltonian
!
                           INTERFACE lanczos_convergence_test
             MODULE PROCEDURE lanczos_convergence_test_d,          &
                              lanczos_convergence_test_z  
                       END INTERFACE lanczos_convergence_test
!
                           INTERFACE preconditioned_lanczos
             MODULE PROCEDURE Preconditioned_lanczos_d,                           &
                              preconditioned_lanczos_z
                       END INTERFACE preconditioned_lanczos
!
                           INTERFACE preconditioned_lanczos_vec
             MODULE PROCEDURE preconditioned_lanczos_vec_d,                       &
                              preconditioned_lanczos_vec_z
                       END INTERFACE preconditioned_lanczos_vec
!
                           INTERFACE preconditioned_lanczos_hamiltonian
             MODULE PROCEDURE preconditioned_lanczos_hamiltonian
                       END INTERFACE preconditioned_lanczos_hamiltonian
!
                           INTERFACE preconditioned_lanczos_convergence_test
             MODULE PROCEDURE preconditioned_lanczos_convergence_test_d,          &
                              preconditioned_lanczos_convergence_test_z  
                       END INTERFACE preconditioned_lanczos_convergence_test
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!**begin prologue     lanczos_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Perform the Lanczos recursion at a given timestep.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       lanc_propagation_d
!***********************************************************************
  SUBROUTINE lanczos_d(v_in,v_out)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                :: v_in
  REAL*8, DIMENSION(:)                :: v_out
  CHARACTER (LEN=5)                   :: itoc
  INTEGER                             :: i
  REAL*8                              :: ddot
!
  null_vec = .false.
  convergence = 'unconverged'
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
  lanczos_tri_d=0
!
!
!           Initial vector
!
  local_time=secnds(0.0)
  anorm = 1.d0 / sqrt (ddot(n3d , v_in(:) , 1, v_in(:) , 1 ))
  vec_d(:,0) = anorm * v_in(:)
  v_in(:) = anorm * vec_d(:,0)
  total_time(1) = secnds(0.0) - local_time + total_time(1)
  DO i=0, maxvec-1
!
     WRITE(iout,3) i


!
!           Begin Lanczos iterations
!
!
!            Lets look at the input vector
!
     IF(log_iterative(10).or.log_iterative(2)) THEN
        title='Lanczos Vector iteration = '//itoc(i)
        call prntfm(title,vec_d(:,i),n3d,1,n3d,maxvec+1,iout)
     END IF
!
     call lanczos_vec(v_out,vec_d,i)
!
!
!           If the number of iterations has reached the size of the
!           original matrix without breakdown, we have to be converged.
!
    IF( convergence == 'maximal') THEN
        WRITE(iout,4) eig(0)
        v_in(:) =v_out(:)
        norm = ddot(n3d,v_in,1,v_in,1)
        total_number_of_iterations = total_number_of_iterations + i + 1
        return
!
!           We are actually converged to the desired tolerance.
!
    ELSE IF(convergence == 'converged') THEN
        WRITE(iout,5) i  
        v_in(:) =v_out(:) 
        norm = ddot(n3d,v_in,1,v_in,1)
        total_number_of_iterations = total_number_of_iterations + i + 1
        return
    END IF    
!
!           We need more iterations
!
    total_time(8) = secnds(0.0) - local_time + total_time(8)
  END DO   
  write(iout,6) maxvec
  stop
1 FORMAT('***********************************************'  &
         '*************************') 
2 FORMAT(/,20X,'Beginning Lanczos Iterations At Zero')
3 FORMAT(/,20x,'Lanczos Iteration = ',i4)
4 FORMAT(/,20x,'Number of Iterations Maximal.  Energy = ', e15.8)
5 FORMAT(/,20x,'Convergence After ',i5, ' Iterations')
6 FORMAT(/,20x,'No Convergence After ',i5, 'Iterations')
END SUBROUTINE lanczos_d
!***********************************************************************
!***********************************************************************
!**begin prologue     lanczos_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Perform the Lanczos recursion at a given timestep.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       lanczos_z
!***********************************************************************
  SUBROUTINE lanczos_z(v_in,v_out)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)            :: v_in, v_out
  INTEGER                             :: i
  CHARACTER (LEN=5)                   :: itoc
  COMPLEX*16                          :: cdotc
!
  null_vec = .false.
  convergence = 'unconverged'
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
  lanczos_tri_z=0
!
!
! 
  anorm = 1.d0 / sqrt (cdotc(n3d , v_in(:) , 1, v_in(:) , 1 ))
  vec_z(:,0) = anorm * v_in(:)
  v_in(:) = anorm * vec_z(:,0)
  DO i=0, maxvec-1
     WRITE(iout,3) i
!
!
!
!            Lets look at the input vector
!
     IF(log_iterative(10).or.log_iterative(2)) THEN
        title='Lanczos Vector iteration = '//itoc(i)
        call prntcm(title,vec_z(:,i),n3d,1,n3d,maxvec+1,iout)
     END IF
!
    call lanczos_vec(v_out,vec_z,i)
    IF( convergence == 'maximal') THEN
        WRITE(iout,4)
        v_in(:) =v_out(:)
        norm = cdotc(n3d,v_in,1,v_in,1)
        return
    ELSE IF(convergence == 'converged') THEN
        WRITE(iout,5) i  
        v_in(:) =v_out(:)     
        norm = cdotc(n3d,v_in,1,v_in,1)
        return
    END IF    
  END DO   
  write(iout,6) maxvec
  stop
1 FORMAT('***********************************************'  &
         '*************************') 
2 FORMAT(/,5X,'Beginning Lanczos Iterations Beginning At Zero')
3 FORMAT(/,10x,'Lanczos Iteration = ',i4)
4 FORMAT(/,20x,'Number of Iterations Maximal')
5 FORMAT(/,20x,'Convergence After ',i5, ' Iterations')
6 FORMAT(/,20x,'No Convergence After ',i5, 'Iterations')
END SUBROUTINE lanczos_z
!***********************************************************************
!***********************************************************************
!**begin prologue     lanczos_vec_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Perform the Lanczos recursion at a given timestep.
!**                   There is a standard and preconditioned version.
!**                                 Standard                   
!***                  b_{i+1} |V_{i+1}> = [ H - a_{i} ] |V_{i}> - b_{i} |V_{i-1}>
!**                                 Preconditioned
!***                  b_{i+1} S |V_{i+1}> = [ H - a_{i} S ] |V_{i}> 
!                                               -   b_{i} S |V_{i-1}>
!                                         T
!                                   S = LL
!                                    -1     -T  -1
!                                   S  = ( L   L   )
!**references
!**routines called    iosys, util and mdutil
!**end prologue       lanczos_vec_d
!***********************************************************************
  SUBROUTINE lanczos_vec_d(v_out,lanczos_basis,it)
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d,0:maxvec)                :: lanczos_basis
  REAL*8, DIMENSION(n3d)                         :: v_out
  INTEGER                                        :: i, it, lower
  REAL*8                                         :: ddot
  CHARACTER(LEN=4)                               :: itoc
!
!            Note that it runs from 0 to maxvec
!            The starting vector is assumed to be normalized.
!
!            Lets compute the value of H on the input vector in order to
!            begin the generation of the next vector.
!
  local_time=secnds(0.0)
  IF(type == 'general') THEN
     IF (packed) THEN 
         call column_packed_matrix_on_vector                              &
                            (lanczos_basis(1:n3d,it),                     &
                             h_vec_d,                                     &
                             mat_var(3)%matrix_diagonal_d,                &
                             mat_var(3)%packed_columns_d,                 &
                             mat_var(3)%non_zero_columns,                 &
                             mat_var(3)%row_index,n3d)
     ELSE
         call matrix_on_vector( hamiltonian_d, lanczos_basis(:,it), h_vec_d )
     END IF
  ELSE IF (type == 'finite_element') THEN
     CALL finite_element_m_v(lanczos_basis(:,it), h_vec_d )
  END IF
  IF(log_iterative(10).or.log_iterative(3)) THEN
     title='H_on_Lanczos_Vector iteration = '//itoc(it)
     call prntfm(title,h_vec_d,n3d,1,n3d,1,iout)
  END IF
  total_time(2) = secnds(0.0) - local_time + total_time(2)
!
!             Compute the value of a(it) as the
!                  a(i) = <V_i|h|V_i> 
!
  local_time=secnds(0.0)
  a(it) =  ddot(n3d , vec_d(:,it) , 1 , h_vec_d(:) , 1 )
  IF(log_iterative(10).or.log_iterative(1)) THEN
     title='Lanczos Alpha iteration = '//itoc(it)
     call prntfmn(title,a(it),1,1,maxvec+1,1,iout,'e')
  END IF
!
!             Compute eigenvalues and eigenvectors of the projected
!             Hamiltonian
!
  call lanczos_hamiltonian(eig,sub_diagonal,eigen_vectors,it)
  total_time(3) = secnds(0.0) - local_time + total_time(3)
!
!             Test convergence
!
  local_time=secnds(0.0)
  call lanczos_convergence_test(eigen_vectors,lanczos_tri_d,it)
  total_time(4) = secnds(0.0) - local_time + total_time(4)
!
!             Branch based on test.
!
  IF(convergence == 'unconverged') THEN
!
!          Compute next vector at it + 1
!
      lanczos_basis(:,it+1) =  h_vec_d(:) - a(it) * lanczos_basis(:,it)
      IF( it > 0 ) THEN
          lanczos_basis(:,it+1) = lanczos_basis(:,it+1)               &
                                          -                           &
                          b(it) * lanczos_basis(:,it-1)
         IF(log_iterative(10)) THEN
             title='before second schmidt'
             call prntfm(title,lanczos_basis,n3d,it+2,n3d,maxvec+1,iout)
         END IF
!
!           Performing a reorthogonalize orthogonalization can help
!           prevent linear dependence.  Its an option if needed.
!
         local_time=secnds(0.0)
         IF(orthogonalize=='double_schmidt') THEN
            lower = it - 1
         ELSE IF(orthogonalize=='full') THEN
             lower = 0
         ELSE IF(orthogonalize=='partial') THEN
             lower = min(it_min,it-1)
         END IF
         DO i=lower,it
            lanczos_basis(:,it+1) = lanczos_basis(:,it+1)           &
                                  -                                 &
            ddot(n3d,lanczos_basis(:,it+1),1,lanczos_basis(:,i),1)  &
                                  *                                 &
                                    lanczos_basis(:,i) 
         END DO
         total_time(5) = secnds(0.0) - local_time + total_time(5)
         IF(log_iterative(10)) THEN
            title='after reorthogonalization'
            call prntfm(title,lanczos_basis,n3d,it+2,n3d,maxvec+1,iout)
         END IF
      END IF
!
!          Compute b(it+1)  b(i) = <V_i|V_i> 
!
      b(it+1) = SQRT( ddot(n3d , lanczos_basis(:,it+1) , 1 ,       &
                                 lanczos_basis(:,it+1) , 1) ) 
      IF(log_iterative(10).or.log_iterative(2)) THEN
         title='Lanczos Beta iteration = '//itoc(it)
         call prntfm(title,b(it+1),1,1,maxvec+1,maxvec+1,iout)
      END IF
!
!          Test for breakdown.  It occurs if b(i) becomes zero.
!
      IF( b(it+1) <= eps ) THEN
          local_time=secnds(0.0)
          null_vec=.true.
          write(iout,1)
          call random_number(lanczos_basis(:,it+1))
          DO i=0,it
             lanczos_basis(:,it+1) = lanczos_basis(:,it+1)         &
                                   -                               &
                                     lanczos_basis(:,i)            &
                                   *                               &
             ddot(n3d , lanczos_basis(:,it+1) , 1 , lanczos_basis(:,i) , 1 )
          END DO
         b(it+1) = SQRT( ddot(n3d , lanczos_basis(:,it+1) , 1 ,    &
                                    lanczos_basis(:,it+1) , 1) ) 
         IF(log_iterative(10).or.log_iterative(2)) THEN
            title='New Lanczos Beta iteration = '//itoc(it+1)
            call prntfm(title,b(it+1),1,1,maxvec+1,maxvec+1,iout)
         END IF
         lanczos_basis(:,it+1)=lanczos_basis(:,it+1)/b(it+1)
         total_time(6) = secnds(0.0) - local_time + total_time(6)
         b(it+1) = 0.d0
      ELSE
!
!          Normalize vector
!            
          lanczos_basis(:,it+1)=lanczos_basis(:,it+1)/b(it+1)
!
!
      END IF
  END IF
1 FORMAT(/,10x,'Lanczos Breakdown.  Generate New Vector and Continue')
END SUBROUTINE lanczos_vec_d
!***********************************************************************
!***********************************************************************
!**begin prologue     lanczos_vec_z
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
!**end prologue       lanczos_vec_z
!***********************************************************************
  SUBROUTINE lanczos_vec_z(v_out,lanczos_basis,it)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d,0:maxvec)            :: lanczos_basis
  COMPLEX*16, DIMENSION(n3d)                     :: v_out
  INTEGER                                        :: i, it, lower
  COMPLEX*16                                     :: cdotc
  CHARACTER(LEN=4)                               :: itoc
  REAL*8                                         :: t_ran
!
!            Note that it runs from 0 to maxvec
!            The starting vector is assumed to be normalized.
!
!            Lets compute the value of H on the input vector in order to
!            calculate the next a.
!
  IF(type == 'general') THEN
     IF (packed) THEN 
         call column_packed_matrix_on_vector                              &
                            (lanczos_basis(1:n3d,it),                     &
                             h_vec_z,                                     &
                             mat_var(3)%matrix_diagonal_z,                &
                             mat_var(3)%packed_columns_z,                 &
                             mat_var(3)%non_zero_columns,                 &
                             mat_var(3)%row_index,n3d)
     ELSE
         call matrix_on_vector( hamiltonian_z, lanczos_basis(:,it), h_vec_z )
     END IF
  ELSE IF (type == 'finite_element') THEN
      CALL finite_element_m_v(lanczos_basis(:,it), h_vec_z )
  END IF
  IF(log_iterative(10).or.log_iterative(3)) THEN
     title='H_on_Lanczos_Vector iteration = '//itoc(it)
     call prntcmn(title,h_vec_z,n3d,1,n3d,1,iout,'e')
  END IF
!
!             Compute the value of a(it) as the
!                  a(i) = <V_i|h|V_i> 
!
  a(it) =  cdotc(n3d , vec_z(:,it) , 1 , h_vec_z(:) , 1 )
  IF(log_iterative(10).or.log_iterative(1)) THEN
     title='Lanczos Alpha iteration = '//itoc(it)
     call prntfmn(title,a(it),1,1,maxvec+1,1,iout,'e')
  END IF
!
!             Compute eigenvalues and eigenvectors of the projected
!             Hamiltonian
!
  call lanczos_hamiltonian(eig,sub_diagonal,eigen_vectors,it)
!
!            Test convergence
!
  call lanczos_convergence_test(eigen_vectors,lanczos_tri_z,it)
!
  IF(convergence == 'unconverged') THEN
!
!          Compute next vector at it + 1
!
      lanczos_basis(:,it+1) =     h_vec_z(:) - a(it) * lanczos_basis(:,it)
      IF( it > 0 ) THEN
          lanczos_basis(:,it+1) = lanczos_basis(:,it+1)               &
                                          -                           &
                          b(it) * lanczos_basis(:,it-1)
         IF(log_iterative(10)) THEN
            title='before second schmidt'
            call prntcmn(title,lanczos_basis,n3d,it+2,n3d,maxvec+1,iout,'e')
         END IF
!
!           Performing a reorthogonalize orthogonalization can help
!           prevent linear dependence.  Its an option if needed.
!
         IF(orthogonalize=='double_schmidt') THEN
            lower = it - 1
         ELSE IF(orthogonalize=='full') THEN
             lower = 0
         ELSE IF(orthogonalize=='partial') THEN
             lower = min(it_min,it-1)
         END IF
         DO i=lower,it
            lanczos_basis(:,it+1) = lanczos_basis(:,it+1)           &
                                  -                                 &
            cdotc(n3d,lanczos_basis(:,it+1),1,lanczos_basis(:,i),1) &
                                  *                                 &
                                              lanczos_basis(:,i) 
         END DO
         IF(log_iterative(10)) THEN
            title='after reorthogonalization'
            call prntfm(title,lanczos_basis,n3d,it+2,n3d,maxvec+1,iout)
         END IF
      END IF
!
!          Compute b(it+1)
!
      b(it+1) = SQRT( cdotc(n3d , lanczos_basis(:,it+1) , 1 ,       &
                              
                                  lanczos_basis(:,it+1) , 1) ) 
      IF(log_iterative(10).or.log_iterative(2)) THEN
         title='Lanczos Beta iteration = '//itoc(it)
         call prntfm(title,b(it+1),1,1,maxvec+1,maxvec+1,iout)
      END IF
!
!          Test for breakdown
!
      IF( b(it+1) <= eps ) THEN
          null_vec=.true.
          write(iout,1)
          DO i=1,n3d
             call random_number(t_ran)
             lanczos_basis(i,it+1) = t_ran
          END DO
          DO i=0,it
             lanczos_basis(:,it+1) = lanczos_basis(:,it+1)         &
                                   -                               &
                                     lanczos_basis(:,i)            &
                                   *                               &
             cdotc(n3d , lanczos_basis(:,it+1) , 1 , lanczos_basis(:,i) , 1 )
          END DO
          b(it+1) = SQRT( cdotc(n3d , lanczos_basis(:,it+1) , 1 , lanczos_basis(:,it+1) , 1) ) 
          IF(log_iterative(10).or.log_iterative(2)) THEN
            title='New Lanczos Beta iteration = '//itoc(it+1)
             call prntfm(title,b(it+1),1,1,maxvec+1,maxvec+1,iout)
          END IF
          lanczos_basis(:,it+1)=lanczos_basis(:,it+1)/b(it+1)
          b(it+1) = 0.d0
      ELSE
!
!          Normalize vector
!            
          lanczos_basis(:,it+1)=lanczos_basis(:,it+1)/b(it+1)
!
!
      END IF
  END IF
1 FORMAT(/,10x,'Lanczos Breakdown.  Generate New Vector and Continue')
END SUBROUTINE lanczos_vec_z
!***********************************************************************
!***********************************************************************
!**begin prologue     lanczos_hamiltonian
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
!**end prologue       lanczos_hamiltonian
!***********************************************************************
!***********************************************************************
  SUBROUTINE lanczos_hamiltonian(eig,sub_diagonal,eigen_vectors,it)
  IMPLICIT NONE
  REAL*8, DIMENSION(0:maxvec)                    :: eig
  REAL*8, DIMENSION(0:maxvec)                    :: sub_diagonal  
  REAL*8, DIMENSION(0:maxvec,0:maxvec)           :: eigen_vectors  
  INTEGER                                        :: it
  CHARACTER(LEN=4)                               :: itoc
!
!               
  eig(0:it) = a(0:it)
  sub_diagonal(0:it-1) = b(1:it)
  IF(log_iterative(10).or.log_iterative(5)) THEN
     title = 'Diagonal Element Lanczos H iteration = '//itoc(it)
     call prntfmn(title,eig,it+1,1,maxvec+1,1,iout,'e')
     title = 'Off-Diagonal Element Lanczos H iteration = '//itoc(it)
     call prntfmn(title,sub_diagonal,it,1,maxvec+1,1,iout,'e')
  END IF
!
!            Lets get the eigenvalues.
!
  call dstev('v',it+1,eig,sub_diagonal,eigen_vectors,maxvec+1,rwork,info)
  title='Lanczos eigenvalues iteration = '//itoc(it)
  call prntfmn(title,eig,it+1,1,maxvec+1,maxvec+1,iout,'e')
  IF(log_iterative(10).or.log_iterative(7)) THEN
     title='Eigenvectors in lanczos basis iteration = '//itoc(it)
     call prntfmn(title,eigen_vectors,it+1,it+1,maxvec+1,maxvec+1,iout,'e')
  END IF
!
END SUBROUTINE lanczos_hamiltonian
!***********************************************************************
!***********************************************************************
!**begin prologue     lanczos_convergence_test_d
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
!**end prologue       lanczos_convergence_test_d
!***********************************************************************
  SUBROUTINE lanczos_convergence_test_d(eigen_vectors,lanczos_tri,it)
  USE dvrprop_global,     local_scratch=> v_scr_d
  IMPLICIT NONE
  REAL*8, DIMENSION(0:maxvec,0:maxvec)           :: eigen_vectors
  REAL*8, DIMENSION(0:maxvec)                    :: lanczos_tri
  INTEGER                                        :: i, it
  REAL*8                                         :: ddot
  CHARACTER(LEN=4)                               :: itoc
!
!               
!             Compute the projection of the time-dependent state on the initial vector,
!             The initial vector is the first Lanczos vector.  So, this just scales
!             the first component of each eigenvector by an exponential.
!
!             C(t) = Sum_{lambda} Exp(-lambda * deltat ) d_{lambda}(0) C_{lambda}
!                                                     +
!             Compute the d_{lambda}(0) = [C_{lambda}]  S C(0)
!
  local_scratch(1:it+1) = eigen_vectors(0,0:it)
  IF(log_iterative(10)) THEN
     title='Projection of Initial Vector on Eigenbasis iteration = '//itoc(it)
     call prntfm(title,local_scratch,it+1,1,n3d,1,iout)
  END IF
! 
!           Here is the scaling.
!
  local_scratch(1:it+1) = EXP(-eig(0:it) * deltat/hbar) * local_scratch(1:it+1)
  IF(log_iterative(10)) THEN
      title='Projection Multiplied by Exponential iteration = '//itoc(it)
      call prntfm(title,local_scratch,it+1,1,maxvec+1,1,iout)
  END IF
!
!  
!             Now express things in the Lanczos basis instead of the eigenfunction basis.
!
  call ebcx(work_d,maxvec+1,eigen_vectors,maxvec+1,local_scratch,maxvec+1,it+1,it+1,1)
!
!             At this point, if we were converged, we would transform to the
!             original basis and return.  However, if we are not converged
!             we do not have to do that.  So, lets test convergence.
!
!
  IF( it+1 == n3d) THEN
!
!             We have maximized the number of iterations, we are totally finished.  
!             Take the lowest eigenvalue and eigenvector, transform to the original basis and return,
!               
       convergence = 'maximal'
       title='Final Value of Lowest Eigenvalue'
       call prntfm(title,eig,1,1,maxvec+1,1,iout)
       eig_old = eig(0)
       call ebcx(local_scratch,n3d,vec_d,n3d,work_d,maxvec+1,n3d,it+1,1)
       title='Final Value of Exp(-H*t) on Initial Vector'
       call prntfm(title,local_scratch,n3d,1,n3d,1,iout)
       call ebcx(local_scratch,n3d,vec_d,n3d,eigen_vectors,maxvec+1,n3d,it+1,1)
       IF(log_iterative(10).or.log_iterative(9)) THEN     
          title='Final Solution'
          call prntfm(title,local_scratch,n3d,1,n3d,1,iout)
       END IF 
       RETURN
  END IF
  IF (it == 0) THEN
      lanczos_tri(it) = work_d(it)
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          title='lanczos_tri iteration = '//itoc(it)
          call prntfm(title,lanczos_tri,it+1,1,maxvec+1,1,iout)
       END IF
  ELSE
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          title='lanczos_tri iteration = '//itoc(it)
          call prntfm(title,lanczos_tri,it+1,1,maxvec+1,1,iout)
          title='work iteration = '//itoc(it)
          call prntfm(title,work_d,it+1,1,maxvec+1,1,iout)
       END IF
       lanczos_tri(0:it) = lanczos_tri(0:it) - work_d(0:it)
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          title='lanczos_tri iteration = '//itoc(it)
          call prntfm(title,lanczos_tri,it+1,1,maxvec+1,1,iout)
       END IF
  END IF
  wfn_tst = sqrt ( ddot( it+1 , lanczos_tri , 1 ,lanczos_tri , 1 ) )
  write(iout,1) wfn_tst
  lanczos_tri(0:it) = work_d(0:it)
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
       title='Converged Value of Lowest Eigenvalue'
       call prntfm(title,eig,1,1,maxvec+1,1,iout)
!
!    Transform to original basis
!    C_alpha(t) = Sum vec_d(alpha,i) * A(i)
!    C(t) = Sum C_alpha(t) * q_alpha
!
     eig_old = eig(0)
     call ebcx(local_scratch,n3d,vec_d,n3d,eigen_vectors,maxvec+1,n3d,it+1,1)
     IF(log_iterative(10).or.log_iterative(9)) THEN     
        title='Final Solution'
        call prntfm(title,local_scratch,n3d,1,n3d,1,iout)
     END IF
     RETURN
  END IF
1 FORMAT(/,10x,'RMS = sqrt(< ( v_{i} - v_{i-1}) |( v_{i} - v_{i-1}) >) = ',e15.8)
END SUBROUTINE lanczos_convergence_test_d
!***********************************************************************
!***********************************************************************
!**begin prologue     lanczos_convergence_test_z
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
!**end prologue       lanczos_convergence_test_z
!***********************************************************************
  SUBROUTINE lanczos_convergence_test_z(eigen_vectors,lanczos_tri,it)
  USE dvrprop_global,     local_scratch=> v_scr_z
  IMPLICIT NONE
  REAL*8, DIMENSION(0:maxvec,0:maxvec)               :: eigen_vectors
  COMPLEX*16, DIMENSION(0:maxvec)                    :: lanczos_tri
  INTEGER                                            :: i, it
  COMPLEX*16                                         :: cdotc
  CHARACTER(LEN=4)                                   :: itoc
!
!               
  local_scratch(1:it+1) = eigen_vectors(0,0:it)
  IF(log_iterative(10)) THEN
     title='Projection of Initial Vector on Eigenbasis iteration = '//itoc(it)
     call prntcmn(title,local_scratch,it+1,1,n3d,1,iout,'e')
  END IF
  local_scratch(1:it+1) = EXP(-eye*eig(0:it) * deltat/hbar) * local_scratch(1:it+1)
  IF(log_iterative(10)) THEN
     title='Projection Multiplied by Exponential iteration = '//itoc(it)
     call prntcmn(title,local_scratch,it+1,1,maxvec+1,1,iout,'e')
  END IF
!
!             Transform to Lanczos basis.
!
  call ebccx(work_z,maxvec+1,eigen_vectors,maxvec+1,local_scratch,n3d,it+1,it+1,1)
!
!             At this point, if we were converged, we would transform to the
!             original basis and return.  However, if we are not converged
!             we do not have to do that.  So, lets test convergence.
  IF( it+1 == n3d) THEN
!
       convergence = 'maximal'
       call cebcx(local_scratch,n3d,vec_z,n3d,work_z,maxvec+1,n3d,it+1,1)
!       local_scratch =local_scratch/anorm
       title='Value of Exp(-i*H*t) on Initial Vector Final'
       call prntcmn(title,local_scratch,n3d,1,n3d,1,iout,"e")
       IF(log_iterative(10).or.log_iterative(9)) THEN     
          title='Final Solution'
          call prntcmn(title,local_scratch,n3d,1,n3d,1,iout,'e')
       END IF 
       RETURN
  END IF
  IF (it == 0) THEN
      lanczos_tri(it) = work_z(it)
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          title='lanczos_tri iteration = '//itoc(it)
          call prntcmn(title,lanczos_tri,it+1,1,maxvec+1,1,iout,'e')
       END IF
  ELSE
      IF(log_iterative(10).or.log_iterative(8)) THEN     
         title='lanczos_tri iteration = '//itoc(it)
         call prntcmn(title,lanczos_tri,it+1,1,maxvec+1,1,iout,'e')
         title='work iteration = '//itoc(it)
         call prntcmn(title,work_z,it+1,1,maxvec+1,1,iout,'e')
      END IF
      lanczos_tri(0:it) = lanczos_tri(0:it) - work_z(0:it)
      IF(log_iterative(10).or.log_iterative(8)) THEN     
         title='lanczos_tri iteration = '//itoc(it)
         call prntcmn(title,lanczos_tri,it+1,1,maxvec+1,1,iout,'e')
      END IF
  END IF
  wfn_tst = sqrt ( cdotc( it+1 , lanczos_tri , 1 ,lanczos_tri , 1 ) )
  write(iout,1) wfn_tst
  lanczos_tri(0:it) = work_z(0:it)
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
     call cebcx(local_scratch,n3d,vec_z,n3d,work_z,maxvec+1,n3d,it+1,1)
!     title='final solution before renormalization'
!     call prntfm(title,local_scratch,n3d,1,n3d,1,iout)
!     local_scratch = local_scratch / anorm
!     title='final solution after renormalization'
!     call prntfm(title,local_scratch,n3d,1,n3d,1,iout)
     IF(log_iterative(10).or.log_iterative(9)) THEN     
        title='Final Solution'
        call prntcm(title,local_scratch,n3d,1,n3d,1,iout,'e')
     END IF
     RETURN
  END IF
1 FORMAT(/,10x,'RMS = sqrt(< ( v_{i} - v_{i-1}) |( v_{i} - v_{i-1}) >) = ',e15.8)
END SUBROUTINE lanczos_convergence_test_z
!***********************************************************************
!***********************************************************************
!**begin prologue     preconditioned_lanczos_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Perform the Lanczos recursion at a given timestep.
!**description        
!**                   
!**references
!**routines called    
!**end prologue       preconditioned_lanczos_d
!***********************************************************************
  SUBROUTINE preconditioned_lanczos_d(v_in,v_out)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                           :: v_in
  REAL*8, DIMENSION(:)                           :: v_out
  CHARACTER (LEN=5)                              :: itoc
  INTEGER                                        :: i
  REAL*8                                         :: ddot
!
  null_vec = .false.
  convergence = 'unconverged'
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
  lanczos_tri_d=0
!
!
!           Compute Initial S-normalized vector
!
  local_time=secnds(0.0)
  IF (packed) THEN 
      call column_packed_matrix_on_vector                              &
                         (v_in(:),                                     &
                          s_vec_d(:,0),                                &
                          mat_var(4)%matrix_diagonal_d,                &
                          mat_var(4)%packed_columns_d,                 &
                          mat_var(4)%non_zero_columns,                 &
                          mat_var(4)%row_index,n3d)
  ELSE
      call ebc(s_vec_d(:,0), overlap_d(:,:), v_in(:), n3d, n3d,1)
  END IF
  anorm = 1.d0 / sqrt (ddot(n3d , v_in(:) , 1, s_vec_d(:,0) , 1 ))
  vec_d(:,0) = anorm * v_in(:)
  s_vec_d(:,0) = anorm * s_vec_d(:,0)
  v_in = anorm * v_in
  total_time(1) = secnds(0.0) - local_time + total_time(1)
  DO i=0, maxvec-1
!
     WRITE(iout,3) i
!
!
!            Lets look at the input vector
!
     IF(log_iterative(10).or.log_iterative(2)) THEN
        title='Lanczos Vector iteration = '//itoc(i)
        call prntfmn(title,vec_d(:,i),n3d,1,n3d,maxvec+1,iout,'e')
     END IF
!
!           Begin Lanczos iterations
!
     call preconditioned_lanczos_vec(v_out,vec_d,s_vec_d,i)
!
!
!           If the number of iterations has reached the size of the
!           original matrix without breakdown, we have to be converged.
!
    local_time = secnds(0.0)
    IF( convergence == 'maximal') THEN
        WRITE(iout,4) eig(0)
        IF (packed) THEN 
            call column_packed_matrix_on_vector                              &
                               (v_out(:),                                    &
                                s_vec_d(:,0),                                &
                                mat_var(4)%matrix_diagonal_d,                &
                                mat_var(4)%packed_columns_d,                 &
                                mat_var(4)%non_zero_columns,                 &
                                mat_var(4)%row_index,n3d)
        ELSE
            call ebc(s_vec_d(:,0), overlap_d(:,:), v_out(:), n3d, n3d,1)
        END IF
        norm = 1.d0 / sqrt (ddot(n3d , v_out(:) , 1, s_vec_d(:,0) , 1 ))
        v_in = v_out
        total_number_of_iterations = total_number_of_iterations + i + 1
        return
!
!           We are actually converged to the desired tolerance.
!
    ELSE IF(convergence == 'converged') THEN
        WRITE(iout,5) i  
        IF (packed) THEN 
            call column_packed_matrix_on_vector                              &
                               (v_out(:),                                    &
                                s_vec_d(:,0),                                &
                                mat_var(4)%matrix_diagonal_d,                &
                                mat_var(4)%packed_columns_d,                 &
                                mat_var(4)%non_zero_columns,                 &
                                mat_var(4)%row_index,n3d)
        ELSE
            call ebc(s_vec_d(:,0), overlap_d(:,:), v_out(:), n3d, n3d,1)
        END IF
        norm = 1.d0 / sqrt (ddot(n3d , v_out(:) , 1, s_vec_d(:,0) , 1 ))
        v_in = v_out
        total_number_of_iterations = total_number_of_iterations + i + 1
        return
    END IF    
    total_time(8) = secnds(0.0) - local_time + total_time(8)
!
!           We need more iterations
!
  END DO   
  write(iout,6) maxvec
  stop
1 FORMAT('***********************************************'  &
         '*************************') 
2 FORMAT(/,20X,'Beginning Lanczos Iterations At Zero')
3 FORMAT(/,20x,'Lanczos Iteration = ',i4)
4 FORMAT(/,20x,'Number of Iterations Maximal.  Energy = ', e15.8)
5 FORMAT(/,20x,'Convergence After ',i5, ' Iterations')
6 FORMAT(/,20x,'No Convergence After ',i5, 'Iterations')
END SUBROUTINE preconditioned_lanczos_d
!***********************************************************************
!***********************************************************************
!**begin prologue     preconditioned_lanczos_z
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Perform the Lanczos recursion at a given timestep.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       preconditioned_lanczos_z
!***********************************************************************
  SUBROUTINE preconditioned_lanczos_z(v_in,v_out)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)            :: v_in
  COMPLEX*16, DIMENSION(:)            :: v_out
  CHARACTER (LEN=5)                   :: itoc
  INTEGER                             :: i
  COMPLEX*16                          :: cdotc
!
  null_vec = .false.
  convergence = 'unconverged'
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
  lanczos_tri_z=0
!
!
! 
  IF (packed) THEN 
      call column_packed_matrix_on_vector                              &
                         (v_in(:),                                     &
                          s_vec_z(:,0),                                &
                          mat_var(4)%matrix_diagonal_z,                &
                          mat_var(4)%packed_columns_z,                 &
                          mat_var(4)%non_zero_columns,                 &
                          mat_var(4)%row_index,n3d)
  ELSE
      call cebc(s_vec_z(:,0), overlap_z(:,:), v_in(:), n3d, n3d,1)
  END IF
  anorm = 1.d0/sqrt(cdotc(n3d,v_in,1,s_vec_z(:,0),1))
  vec_z(:,0) = anorm * v_in(:)
  s_vec_z(:,0) = anorm * s_vec_z(:,0)
  v_in = anorm * v_in
!
  DO i=0, maxvec-1
     WRITE(iout,3) i
!
!
!
!            Lets look at the input vector
!
     IF(log_iterative(10).or.log_iterative(2)) THEN
        title='Lanczos Vector iteration = '//itoc(i)
        call prntcmn(title,vec_z(:,i),n3d,1,n3d,maxvec+1,iout,'e')
     END IF
!
!           Begin Lanczos iterations
!
    call preconditioned_lanczos_vec(v_out,vec_z,s_vec_z,i)
!
!
!           If the number of iterations has reached the size of the
!           original matrix without breakdown, we have to be converged.
!
    IF( convergence == 'maximal') THEN
        WRITE(iout,4)
        v_in(:) =v_out(:)
        IF (packed) THEN 
            call column_packed_matrix_on_vector                              &
                               (v_in(:),                                     &
                                s_vec_z(:,0),                                &
                                mat_var(4)%matrix_diagonal_z,                &
                                mat_var(4)%packed_columns_z,                 &
                                mat_var(4)%non_zero_columns,                 &
                                mat_var(4)%row_index,n3d)
        ELSE
            call cebc(s_vec_z(:,0), overlap_z(:,:), v_in(:), n3d, n3d,1)
        END IF
        norm = sqrt (cdotc(n3d , v_in(:) , 1, s_vec_z(:,0) , 1 ))
        IF(compute_energy) THEN 
           IF(type == 'general') THEN
              IF (packed) THEN 
                  call column_packed_matrix_on_vector                        &
                               (v_in(:),                                     &
                                h_vec_z(:),                                  &
                                mat_var(3)%matrix_diagonal_z,                &
                                mat_var(3)%packed_columns_z,                 &
                                mat_var(3)%non_zero_columns,                 &
                                mat_var(3)%row_index,n3d)
              ELSE
                  call matrix_on_vector( hamiltonian_z, v_in, h_vec_z )
              END IF
           ELSE IF (type == 'finite_element') THEN
              CALL finite_element_m_v(v_in, h_vec_z )
           END IF
           eig(1) = cdotc(n3d,v_in,1,h_vec_z,1)
           write(iout,*) 'Expectation Value of Hamiltonian = ', eig(1)
        END IF
        return
    ELSE IF(convergence == 'converged') THEN
        WRITE(iout,5) i  
        v_in(:) = v_out(:)     
        IF (packed) THEN 
            call column_packed_matrix_on_vector                              &
                               (v_in(:),                                     &
                                s_vec_z(:,0),                                &
                                mat_var(4)%matrix_diagonal_z,                &
                                mat_var(4)%packed_columns_z,                 &
                                mat_var(4)%non_zero_columns,                 &
                                mat_var(4)%row_index,n3d)
        ELSE
            call cebc(s_vec_z(:,0), overlap_z(:,:), v_in(:), n3d, n3d,1)
        END IF
        norm = sqrt (cdotc(n3d , v_in(:) , 1, s_vec_z(:,0) , 1 ))
        IF(compute_energy) THEN 
           IF(type == 'general') THEN
              IF (packed) THEN 
                  call column_packed_matrix_on_vector                        &
                               (v_in(:),                                     &
                                h_vec_z(:),                                  &
                                mat_var(3)%matrix_diagonal_z,                &
                                mat_var(3)%packed_columns_z,                 &
                                mat_var(3)%non_zero_columns,                 &
                                mat_var(3)%row_index,n3d)
              ELSE
                  call matrix_on_vector( hamiltonian_z, v_in, h_vec_z )
              END IF
           ELSE IF (type == 'finite_element') THEN
              CALL finite_element_m_v(v_in, h_vec_z )
           END IF
           eig(1) = cdotc(n3d,v_in,1,h_vec_z,1)
           write(iout,*) 'Expectation Value of Hamiltonian = ', eig(1)
        END IF
        return
    END IF    
  END DO   
  write(iout,6) maxvec
  stop
1 FORMAT('***********************************************'  &
         '*************************') 
2 FORMAT(/,5X,'Beginning Lanczos Iterations Beginning At Zero')
3 FORMAT(/,10x,'Lanczos Iteration = ',i4)
4 FORMAT(/,20x,'Number of Iterations Maximal')
5 FORMAT(/,20x,'Convergence After ',i5, ' Iterations')
6 FORMAT(/,20x,'No Convergence After ',i5, 'Iterations')
END SUBROUTINE preconditioned_lanczos_z
!***********************************************************************
!***********************************************************************
!**begin prologue     preconditioned_lanczos_vec_d
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Perform the Lanczos recursion at a given timestep.
!***                  b_{i+1} S |V_{i+1}> = [ H - a_{i} S ] |V_{i}> 
!                                               -   b_{i} S |V_{i-1}>
!                                        T
!                                   S = U U
!                                    -1     -1  -T
!                                   S  = ( U   U   )
!**description        The Lanczos recursion for a preconditioned matrix requires
!**                   not only the vectors but the value of the preconditioner
!**                   on these vectors and the solution os some linear systems.
!**                   Actually all that needs to be stored in central memory are
!**                   the values of these for three iterations.  Since we need all
!**                   of the vectors later, we keep these in memory but only store
!**                   the previous two iterations for the operation of the preconditioner
!**                   on the vectors.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       preconditioned_lanczos_vec_d
!***********************************************************************
  SUBROUTINE preconditioned_lanczos_vec_d(v_out,lanczos_basis,         &
                                          s_on_lanczos_basis,it)
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d,0:maxvec)                :: lanczos_basis
  REAL*8, DIMENSION(n3d,0:maxvec)                :: s_on_lanczos_basis
  REAL*8, DIMENSION(n3d)                         :: v_out
  INTEGER                                        :: i
  INTEGER                                        :: it
  INTEGER                                        :: low
  REAL*8                                         :: ddot
  CHARACTER(LEN=4)                               :: itoc
  REAL*8                                         :: over
!
!            Note that it runs from 0 to maxvec
!            The starting vector is assumed to be S normalized.
!
!
!            Lets compute the value of H on the input vector in order to
!            generate the next vector.
!
  local_time=secnds(0.0)
  IF(type == 'general') THEN
     IF (packed) THEN 
         call column_packed_matrix_on_vector                              &
                            (lanczos_basis(1:n3d,it),                     &
                             h_vec_d,                                     &
                             mat_var(3)%matrix_diagonal_d,                &
                             mat_var(3)%packed_columns_d,                 &
                             mat_var(3)%non_zero_columns,                 &
                             mat_var(3)%row_index,n3d)
     ELSE
         call matrix_on_vector( hamiltonian_d, lanczos_basis(:,it), h_vec_d )
     END IF
  ELSE IF (type == 'finite_element') THEN
     CALL finite_element_m_v(lanczos_basis(:,it), h_vec_d )
  END IF
  total_time(2) = secnds(0.0) - local_time + total_time(2)
  IF(log_iterative(10).or.log_iterative(3)) THEN
     title='H_on_Lanczos_Vector iteration = '//itoc(it)
     call prntfmn(title,h_vec_d,n3d,1,n3d,1,iout,'e')
  END IF
!
!             Compute the value of a(it) as the
!                  a(i) = <V_i|h|V_i> 
!
  local_time=secnds(0.0)
  a(it) =  ddot(n3d , vec_d(:,it) , 1 , h_vec_d(:) , 1 )
  IF(log_iterative(10).or.log_iterative(1)) THEN
     title='Lanczos Alpha iteration = '//itoc(it)
     call prntfmn(title,a(it),1,1,maxvec+1,1,iout,'e')
  END IF
!
!             Compute eigenvalues and eigenvectors of the projected
!             Hamiltonian
!
  call preconditioned_lanczos_hamiltonian(eig,sub_diagonal,eigen_vectors,it)
  total_time(3) = secnds(0.0) - local_time + total_time(3)
!
!             Test convergence
!
  local_time=secnds(0.0)
  call preconditioned_lanczos_convergence_test(eigen_vectors,lanczos_tri_d,it)
  total_time(4) = secnds(0.0) - local_time + total_time(4)
!
!             Branch based on test.
!
  IF(convergence == 'unconverged') THEN
!
!          Compute next vector at it + 1
!
!                 First compute
!
!    v_out =  b_{i+1} S |V_{i+1}> = [ H - a_{i} S ] |V_{i}> 
!                                               -   b_{i} S |V_{i-1}>
     v_out(:) = h_vec_d(:)  - a(it) * s_on_lanczos_basis(:,it)
     IF( it > 0 ) THEN
!
!        v_out is proportional to S|V(it+1)>
!
         v_out(:) =v_out(:) - b(it) * s_on_lanczos_basis(:,it-1)
!
         IF(log_iterative(10)) THEN
            title='before reorthogonalization'
            call prntfmn(title,v_out,n3d,1,n3d,1,iout,'e')
         END IF
!
!        Perform an S reorthogonalization to  prevent linear dependence.
!
         local_time=secnds(0.0)
         IF(orthogonalize=='double_schmidt') THEN
             low = it - 1
         ELSE IF(orthogonalize=='full') THEN
             low = 0
         ELSE IF(orthogonalize=='partial') THEN
            low = min(it_min,it-1)
         END IF
         DO i=low, it
            v_out(:) = v_out(:)                                              &
                                 -                                           &
                       ddot(n3d,v_out(:),1,lanczos_basis(:,i),1 )            & 
                                 * s_on_lanczos_basis(:,i)
         END DO
!
         total_time(5) = secnds(0.0) - local_time + total_time(5)
         IF(log_iterative(10)) THEN
            title='after reorthogonalization'
            call prntfmn(title,v_out,n3d,1,n3d,1,iout,'e')
         END IF
     END IF
!
!        To compute the next beta, we use the formula,
!
!        2           -1 -T             -T        -T
!    beta  = <v_out| U U   |v_out) = < U v_out | U v_out >
!
!        This comes from the S normalization requirement on the new vector.
!                                  -T
!        Compute the solution to  U  |x> = |v_out>
!        in order to calculate the next beta.
!
     local_time=secnds(0.0)
     IF (packed) THEN
         Call Packed_Solve (lanczos_basis(:,it+1),                           &
                            v_out,                                           &
                            mat_var(2)%non_zero_columns,                     &
                            mat_var(2)%row_index,                            &
                            mat_var(2)%packed_columns_d,                     &
                            mat_var(1)%matrix_diagonal_d,                    &
                            mat_var(2)%number,'t',it,.false.)
     ELSE
         Call Triangular_Solve(lanczos_basis(:,it+1),v_out,'t',it,.false.)
     END IF
!
!
!         Compute b(it+1)  b(i) = Sqrt ( <V_i|V_i> )
!
     b(it+1) = SQRT( ddot(n3d , lanczos_basis(:,it+1) , 1 ,                  &
                                lanczos_basis(:,it+1) , 1) ) 
     IF(log_iterative(10).or.log_iterative(2)) THEN
        title='Lanczos Beta iteration = '//itoc(it)
        call prntfmn(title,b(it+1),1,1,maxvec+1,maxvec+1,iout,'e')
     END IF
!
!          Store the value of the overlap on the new basis function.
!
     s_on_lanczos_basis(:,it+1) = v_out(:) / b(it+1)
!
!          Solve for final new basis function
!
     IF (packed) THEN
         Call Packed_Solve (lanczos_basis(:,it+1),                           &
                            lanczos_basis(:,it+1),                           &
                            mat_var(1)%non_zero_columns,                     &
                            mat_var(1)%row_index,                            &
                            mat_var(1)%packed_columns_d,                     &
                            mat_var(1)%matrix_diagonal_d,                    &
                            mat_var(1)%number,'n',it,.false.)
     ELSE
         Call Triangular_Solve(lanczos_basis(:,it+1),                        &
                               lanczos_basis(:,it+1),                        &
                               'n',it,.false.)
     END IF
     total_time(7) = secnds(0.0) - local_time + total_time(7)
     lanczos_basis(:,it+1) =  lanczos_basis(:,it+1) / b(it+1)
!
!
!        Test for breakdown.  It occurs if b(i) becomes zero.
!
     IF( b(it+1) <= eps ) THEN
         local_time = secnds(0.0)
         null_vec=.true.
!
!        Compute a random vector.
!
         write(iout,1)
         call random_number(lanczos_basis(:,it+1))
         IF (packed) THEN 
             call column_packed_matrix_on_vector                              &
                                (lanczos_basis(:,it+1),                       &
                                 v_out,                                       &
                                 mat_var(4)%matrix_diagonal_d,                &
                                 mat_var(4)%packed_columns_d,                 &
                                 mat_var(4)%non_zero_columns,                 &
                                 mat_var(4)%row_index,n3d)
         ELSE
             call ebc(v_out,overlap_d,lanczos_basis(:,it+1), n3d, n3d, 1)        
         END IF
!
!        S- orthogonalize it to the previous vectors.
!
         DO i=0,it
            v_out(:) = v_out(:)                                              &
                                 -                                           &
                       ddot(n3d,v_out,1,lanczos_basis(:,i),1 )               &
                                 * s_on_lanczos_basis(:,i)
         END DO
!
!        Normalize it.
!
         IF (packed) THEN 
             call column_packed_matrix_on_vector                              &
                                (v_out,                                       &
                                 s_on_lanczos_basis(:,it+1),                  &
                                 mat_var(4)%matrix_diagonal_d,                &
                                 mat_var(4)%packed_columns_d,                 &
                                 mat_var(4)%non_zero_columns,                 &
                                 mat_var(4)%row_index,n3d)
         ELSE
             call ebc(s_on_lanczos_basis(:,it+1), overlap_d, v_out, n3d, n3d, 1)        
         END IF
         b(it+1) = SQRT( ddot(n3d , v_out , 1 , s_on_lanczos_basis(:,it+1) , 1) ) 
         s_on_lanczos_basis(:,it+1) = v_out(:) / ( b(it+1) * b(it+1) )
         IF(log_iterative(10).or.log_iterative(2)) THEN
            title='New Lanczos Beta iteration = '//itoc(it+1)
            call prntfmn(title,b(it+1),1,1,maxvec+1,maxvec+1,iout,'e')
         END IF
         IF (packed) THEN
             Call Packed_Solve (lanczos_basis(:,it+1),                       &
                                s_on_lanczos_basis(:,it+1),                  &
                                mat_var(2)%non_zero_columns,                 &
                                mat_var(2)%row_index,                        &
                                mat_var(2)%packed_columns_d,                 &
                                mat_var(1)%matrix_diagonal_d,                &
                                mat_var(2)%number,'t',it,.false.)
             Call Packed_Solve (lanczos_basis(:,it+1),                       &
                                lanczos_basis(:,it+1),                       &
                                mat_var(1)%non_zero_columns,                 &
                                mat_var(1)%row_index,                        &
                                mat_var(1)%packed_columns_d,                 &
                                mat_var(1)%matrix_diagonal_d,                &
                                mat_var(1)%number,'n',it,.false.)
         ELSE
             Call Triangular_Solve(lanczos_basis(:,it+1),                    &
                                   s_on_lanczos_basis(:,it+1),'t',it,.false.)
             Call Triangular_Solve(lanczos_basis(:,it+1),                    &
                                   lanczos_basis(:,it+1),'n',it,.false.)
         END IF
         total_time(6) = secnds(0.0) - local_time + total_time(6)
         b(it+1) = 0.d0
     END IF
!
!
  END IF
1 FORMAT(/,10x,'Lanczos Breakdown.  Generate New Vector and Continue')
END SUBROUTINE preconditioned_lanczos_vec_d
!***********************************************************************
!***********************************************************************
!**begin prologue     preconditioned_lanczos_vec_z
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
!**end prologue       preconditioned_lanczos_vec_z
!***********************************************************************
  SUBROUTINE preconditioned_lanczos_vec_z(v_out,lanczos_basis,         &
                                          s_on_lanczos_basis,it)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d,0:maxvec)            :: lanczos_basis
  COMPLEX*16, DIMENSION(n3d,0:1)                 :: s_on_lanczos_basis
  COMPLEX*16, DIMENSION(n3d)                     :: v_out
  INTEGER                                        :: i, j
  INTEGER                                        :: it
  INTEGER                                        :: low
  COMPLEX*16                                     :: cdotc
  CHARACTER(LEN=4)                               :: itoc
  REAL*8                                         :: t_ran
!
!            Note that it runs from 0 to maxvec
!            The starting vector is assumed to be normalized.
!
!            Lets compute the value of H on the input vector in order to
!            calculate the next a.
!
  IF(type == 'general') THEN
     IF (packed) THEN 
         call column_packed_matrix_on_vector                              &
                            (lanczos_basis(1:n3d,it),                     &
                             h_vec_z,                                     &
                             mat_var(3)%matrix_diagonal_z,                &
                             mat_var(3)%packed_columns_z,                 &
                             mat_var(3)%non_zero_columns,                 &
                             mat_var(3)%row_index,n3d)
     ELSE
         call matrix_on_vector( hamiltonian_z, lanczos_basis(:,it), h_vec_z )
     END IF
  ELSE IF (type == 'finite_element') THEN
      CALL finite_element_m_v(lanczos_basis(:,it), h_vec_z )
  END IF
  IF(log_iterative(10).or.log_iterative(3)) THEN
     title='H_on_Lanczos_Vector iteration = '//itoc(it)
     call prntcmn(title,h_vec_z,n3d,1,n3d,1,iout,'e')
  END IF
!
!             Compute the value of a(it) as the
!                  a(i) = <V_i|h|V_i> 
!
  a(it) =  cdotc(n3d , vec_z(:,it) , 1 , h_vec_z(:) , 1 )
  IF(log_iterative(10).or.log_iterative(1)) THEN
     title='Lanczos Alpha iteration = '//itoc(it)
     call prntfmn(title,a(it),1,1,maxvec+1,1,iout,'e')
  END IF
!
!             Compute eigenvalues and eigenvectors of the projected
!             Hamiltonian
!
  call preconditioned_lanczos_hamiltonian(eig,sub_diagonal,eigen_vectors,it)
!
!            Test convergence
!
  call preconditioned_lanczos_convergence_test(eigen_vectors,lanczos_tri_z,it)
!
!
  IF(convergence == 'unconverged') THEN
!
!          Compute next vector at it + 1
!
!                 First compute
!
!    v_out =  b_{i+1} S |V_{i+1}> = [ H - a_{i} S ] |V_{i}> 
!                                               -   b_{i} S |V_{i-1}>
!
     v_out(:) = h_vec_z(:)  - a(it) * s_on_lanczos_basis(:,it)
     IF( it > 0 ) THEN
!
!        v_out is proportional to S|V(it+1)>
!
         v_out(:) =v_out(:) - b(it) * s_on_lanczos_basis(:,it-1)
!
         IF(log_iterative(10)) THEN
            title='before reorthogonalization'
            call prntcmn(title,v_out,n3d,1,n3d,1,iout,'e')
         END IF
!
!        Perform an S reorthogonalization to  prevent linear dependence.
!
         IF(orthogonalize=='double_schmidt') THEN
             low = it - 1
         ELSE IF(orthogonalize=='full') THEN
             low = 0
         ELSE IF(orthogonalize=='partial') THEN
            low = min(it_min,it-1)
         END IF
         DO i=low, it
            v_out(:) = v_out(:)                                        &
                                 -                                     &
                       cdotc(n3d,v_out(:),1,lanczos_basis(:,i),1 )     &
                                 * s_on_lanczos_basis(:,i)
         END DO
!
         IF(log_iterative(10)) THEN
            title='after reorthogonalization'
            call prntcmn(title,v_out,n3d,1,n3d,1,iout,'e')
         END IF
     END IF
!
!        To compute the next beta, we use the formula,
!
!        2           -1 -T             -T        -T
!    beta  = <v_out| U U  |v_out) = < U  v_out | U v_out >
!        This comes from the S normalization requirement on the new vector.
!                                  -T
!        Compute the solution to  U  |x> = |v_out>
!        in order to calculate the next beta.
!
     IF (packed) THEN
         Call Packed_Solve (lanczos_basis(:,it+1),                           &
                            v_out,                                           &
                            mat_var(2)%non_zero_columns,                     &
                            mat_var(2)%row_index,                            &
                            mat_var(2)%packed_columns_z,                     &
                            mat_var(1)%matrix_diagonal_z,                    &
                            mat_var(2)%number,'c',it,.false.)
     ELSE
         Call Triangular_Solve(lanczos_basis(:,it+1),v_out,'n',it,.false.)
     END IF
!
!
!         Compute b(it+1)  b(i) = Sqrt ( <V_i|V_i> )
!
     b(it+1) = SQRT( cdotc(n3d , lanczos_basis(:,it+1) , 1 ,       &
                                 lanczos_basis(:,it+1) , 1) ) 
     IF(log_iterative(10).or.log_iterative(2)) THEN
        title='Lanczos Beta iteration = '//itoc(it)
        call prntfmn(title,b(it+1),1,1,maxvec+1,maxvec+1,iout,'e')
     END IF
!
!          Store the value of the overlap on the new basis function.
!
     s_on_lanczos_basis(:,it+1) = v_out(:) / b(it+1)
!
!        Solve for final new basis function
!
     IF (packed) THEN
         Call Packed_Solve (lanczos_basis(:,it+1),                            &
                            lanczos_basis(:,it+1),                            &
                            mat_var(1)%non_zero_columns,                      &
                            mat_var(1)%row_index,                             &
                            mat_var(1)%packed_columns_z,                      &                                    mat_var(1)%matrix_diagonal_z,                     &
                            mat_var(1)%number,'n',it,.false.)
     ELSE
         Call Triangular_Solve(lanczos_basis(:,it+1),lanczos_basis(:,it+1),  &
                               'n',it,.false.)
     END IF
     lanczos_basis(:,it+1) =  lanczos_basis(:,it+1) / b(it+1)
!
!
!        Test for breakdown.  It occurs if b(i) becomes zero.
!
     IF( b(it+1) <= eps ) THEN
         null_vec=.true.
         write(iout,1)
         DO i=1,n3d
            call random_number(t_ran)
            lanczos_basis(i,it+1) = t_ran
         END DO
         IF (packed) THEN 
             call column_packed_matrix_on_vector                              &
                                (lanczos_basis(:,it+1),                       &
                                 v_out,                                       &
                                 mat_var(4)%matrix_diagonal_z,                &
                                 mat_var(4)%packed_columns_z,                 &
                                 mat_var(4)%non_zero_columns,                 &
                                 mat_var(4)%row_index,n3d)
         ELSE
             call cebc(v_out,overlap_z,lanczos_basis(:,it+1), n3d, n3d, 1)        
         END IF
         DO i=0,it
            v_out(:) = v_out(:)                                        &
                                 -                                     &
                       cdotc(n3d,v_out,1,lanczos_basis(:,i),1 )        &
                                 * s_on_lanczos_basis(:,i)
         END DO
         IF (packed) THEN 
             call column_packed_matrix_on_vector                              &
                                (v_out,                                       &
                                 s_on_lanczos_basis(:,it+1),                  &
                                 mat_var(4)%matrix_diagonal_z,                &
                                 mat_var(4)%packed_columns_z,                 &
                                 mat_var(4)%non_zero_columns,                 &
                                 mat_var(4)%row_index,n3d)
         ELSE
             call cebc(s_on_lanczos_basis(:,it+1),overlap_z,v_out,n3d, n3d,1)        
         END IF
         b(it+1) = SQRT( cdotc(n3d , v_out , 1 , s_on_lanczos_basis(:,it+1) , 1) ) 
         s_on_lanczos_basis(:,it+1) = v_out(:) / ( b(it+1) * b(it+1) )
         IF(log_iterative(10).or.log_iterative(2)) THEN
            title='New Lanczos Beta iteration = '//itoc(it+1)
            call prntfmn(title,b(it+1),1,1,maxvec+1,maxvec+1,iout,'e')
         END IF
         IF (packed) THEN
             Call Packed_Solve (lanczos_basis(:,it+1),                     &
                                s_on_lanczos_basis(:,it+1),                &
                                mat_var(2)%non_zero_columns,               &
                                mat_var(2)%row_index,                      &
                                mat_var(2)%packed_columns_z,               &
                                mat_var(1)%matrix_diagonal_z,              &
                                mat_var(2)%number,'c',it,.false.)
             Call Packed_Solve (lanczos_basis(:,it+1),                     &
!                                s_on_lanczos_basis(:,it+1),                &
                                lanczos_basis(:,it+1),                &
                                mat_var(1)%non_zero_columns,               &
                                mat_var(1)%row_index,                      &
                                mat_var(1)%packed_columns_z,               &
                                mat_var(1)%matrix_diagonal_z,              &
                                mat_var(1)%number,'n',it,.false.)
         ELSE
             Call Triangular_Solve(lanczos_basis(:,it+1),                    &
                                   s_on_lanczos_basis(:,it+1),'c',it,.false.)
             Call Triangular_Solve(lanczos_basis(:,it+1),                    &
                                   lanczos_basis(:,it+1),'n',it,.false.)
         END IF
         b(it+1) = 0.d0
     END IF
!
!
  END IF
1 FORMAT(/,10x,'Lanczos Breakdown.  Generate New Vector and Continue')
END SUBROUTINE preconditioned_lanczos_vec_z
!***********************************************************************
!***********************************************************************
!**begin prologue     preconditioned_lanczos_hamiltonian
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
!**end prologue       preconditioned_lanczos_hamiltonian
!***********************************************************************
  SUBROUTINE preconditioned_lanczos_hamiltonian(eig,sub_diagonal,eigen_vectors,it)
  IMPLICIT NONE
  REAL*8, DIMENSION(0:maxvec)                    :: eig
  REAL*8, DIMENSION(0:maxvec)                    :: sub_diagonal  
  REAL*8, DIMENSION(0:maxvec,0:maxvec)           :: eigen_vectors  
  INTEGER                                        :: it
  CHARACTER(LEN=4)                               :: itoc
!
!               
  eig(0:it) = a(0:it)
  sub_diagonal(0:it-1) = b(1:it)
  IF(log_iterative(10).or.log_iterative(5)) THEN
     title = 'Diagonal Element Lanczos H iteration = '//itoc(it)
     call prntfmn(title,eig,it+1,1,maxvec+1,1,iout,'e')
     title = 'Off-Diagonal Element Lanczos H iteration = '//itoc(it)
     call prntfmn(title,sub_diagonal,it,1,maxvec+1,1,iout,'e')
  END IF
!
!            Lets get the eigenvalues.
!
  call dstev('v',it+1,eig,sub_diagonal,eigen_vectors,maxvec+1,rwork,info)
  title='Lanczos eigenvalues iteration = '//itoc(it)
  call prntfmn(title,eig,it+1,1,maxvec+1,maxvec+1,iout,'e')
  IF(log_iterative(10).or.log_iterative(7)) THEN
     title='Eigenvectors in lanczos basis iteration = '//itoc(it)
     call prntfmn(title,eigen_vectors,it+1,it+1,maxvec+1,maxvec+1,iout,'e')
  END IF
!
END SUBROUTINE preconditioned_lanczos_hamiltonian
!***********************************************************************
!***********************************************************************
!**begin prologue     preconditioned_lanczos_convergence_test_d
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
!**end prologue       preconditioned_lanczos_convergence_test_d
!***********************************************************************
  SUBROUTINE preconditioned_lanczos_convergence_test_d(eigen_vectors,lanczos_tri,it)
  USE dvrprop_global,           local_scratch => v_scr_d
  IMPLICIT NONE
  REAL*8, DIMENSION(0:maxvec,0:maxvec)           :: eigen_vectors
  REAL*8, DIMENSION(0:maxvec)                    :: lanczos_tri
  INTEGER                                        :: i, it
  REAL*8                                         :: ddot
  CHARACTER(LEN=4)                               :: itoc
!
!               
!             Compute the projection of the time-dependent state on the initial vector,
!             The initial vector is the first Lanczos vector.  So, this just scales
!             the first component of each eigenvector by an exponential.
!
!             C(t) = Sum_{lambda} Exp(-lambda * deltat ) d_{lambda}(0) C_{lambda}
!                                                     +
!             Compute the d_{lambda}(0) = [C_{lambda}]  S C(0)
!
  local_scratch(1:it+1) = eigen_vectors(0,0:it)
  IF(log_iterative(10)) THEN
     title='Projection of Initial Vector on Eigenbasis iteration = '//itoc(it)
     call prntfmn(title,local_scratch,it+1,1,n3d,1,iout,'e')
  END IF
! 
!           Here is the scaling.
!
  local_scratch(1:it+1) = EXP(-eig(0:it) * deltat/hbar) * local_scratch(1:it+1)
  IF(log_iterative(10)) THEN
      title='Projection Multiplied by Exponential iteration = '//itoc(it)
      call prntfmn(title,local_scratch,it+1,1,maxvec+1,1,iout,'e')
  END IF
!
!  
!             Now express things in the Lanczos basis instead of the eigenfunction basis.
!
  call ebcx(work_d,maxvec+1,eigen_vectors,maxvec+1,local_scratch,maxvec+1,it+1,it+1,1)
!
!             At this point, if we were converged, we would transform to the
!             original basis and return.  However, if we are not converged
!             we do not have to do that.  So, lets test convergence.
!
!
  IF( it+1 == n3d) THEN
!
!             We have maximized the number of iterations, we are totally finished.  
!             Take the lowest eigenvalue and eigenvector, transform to the original basis and return,
!               
       convergence = 'maximal'
       title='Final Value of Lowest Eigenvalue'
       call prntfmn(title,eig,1,1,maxvec+1,1,iout,'e')
       eig_old = eig(0)
       call ebcx(local_scratch,n3d,vec_d,n3d,work_d,maxvec+1,n3d,it+1,1)
       title='Final Value of Exp(-H*t) on Initial Vector'
       call prntfmn(title,local_scratch,n3d,1,n3d,1,iout,'e')
       call ebcx(local_scratch,n3d,vec_d,n3d,eigen_vectors,maxvec+1,n3d,it+1,1)
       IF(log_iterative(10).or.log_iterative(9)) THEN     
          title='Final Solution'
          call prntfmn(title,local_scratch,n3d,1,n3d,1,iout,'e')
       END IF 
       RETURN
  END IF
  IF (it == 0) THEN
       lanczos_tri(it) = work_d(it)
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          title='lanczos_tri iteration = '//itoc(it)
          call prntfmn(title,lanczos_tri,it+1,1,maxvec+1,1,iout,'e')
       END IF
  ELSE
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          title='lanczos_tri iteration = '//itoc(it)
          call prntfmn(title,lanczos_tri,it+1,1,maxvec+1,1,iout,'e')
          title='work iteration = '//itoc(it)
          call prntfmn(title,work_d,it+1,1,maxvec+1,1,iout,'e')
       END IF
       lanczos_tri(0:it) = lanczos_tri(0:it) - work_d(0:it)
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          title='lanczos_tri iteration = '//itoc(it)
          call prntfmn(title,lanczos_tri,it+1,1,maxvec+1,1,iout,'e')
       END IF
  END IF
  wfn_tst = sqrt ( ddot( it+1 , lanczos_tri , 1 ,lanczos_tri , 1 ) )
  write(iout,1) wfn_tst
  lanczos_tri(0:it) = work_d(0:it)
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
       title='Converged Value of Lowest Eigenvalue'
       call prntfmn(title,eig,1,1,maxvec+1,1,iout,'e')
!
!    Transform to original basis
!    C_alpha(t) = Sum vec_d(alpha,i) * A(i)
!    C(t) = Sum C_alpha(t) * q_alpha
!
     eig_old = eig(0)
     call ebcx(local_scratch,n3d,vec_d,n3d,eigen_vectors,maxvec+1,n3d,it+1,1)
     IF(log_iterative(10).or.log_iterative(9)) THEN     
         title='Final Solution'
         call prntfmn(title,local_scratch,n3d,1,n3d,1,iout,'e')
     END IF
     RETURN
  END IF
1 FORMAT(/,10x,'RMS = sqrt(< ( v_{i} - v_{i-1}) |( v_{i} - v_{i-1}) >) = ',e15.8)
END SUBROUTINE preconditioned_lanczos_convergence_test_d
!***********************************************************************
!***********************************************************************
!**begin prologue     preconditioned_lanczos_convergence_test_z
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
!**end prologue       preconditioned_lanczos_convergence_test_z
!***********************************************************************
  SUBROUTINE preconditioned_lanczos_convergence_test_z(eigen_vectors,lanczos_tri,it)
  USE dvrprop_global,     local_scratch=> v_scr_z
  IMPLICIT NONE
  REAL*8, DIMENSION(0:maxvec,0:maxvec)               :: eigen_vectors
  COMPLEX*16, DIMENSION(0:maxvec)                    :: lanczos_tri
  COMPLEX*16                                         :: temp
  INTEGER                                            :: i, it
  COMPLEX*16                                         :: cdotc
  CHARACTER(LEN=4)                                   :: itoc
!
!               
  local_scratch(1:it+1) = eigen_vectors(0,0:it)
  IF(log_iterative(10)) THEN
     title='Projection of Initial Vector on Eigenbasis iteration = '//itoc(it)
     call prntcmn(title,local_scratch,it+1,1,n3d,1,iout,'e')
  END IF
  local_scratch(1:it+1) = EXP(-eye*eig(0:it) * deltat/hbar) * local_scratch(1:it+1)
  IF(log_iterative(10)) THEN
     title='Projection Multiplied by Exponential iteration = '//itoc(it)
     call prntcmn(title,local_scratch,it+1,1,maxvec+1,1,iout,'e')
  END IF
!
!             Transform to Lanczos basis.
!
  call ebccx(work_z,maxvec+1,eigen_vectors,maxvec+1,local_scratch,n3d,it+1,it+1,1)
!
!          At this point, if we were converged, we would transform to the
!          original basis and return.  However, if we are not converged
!          we do not have to do that.  So, lets test convergence.
  IF( it+1 == n3d) THEN
!
       convergence = 'maximal'
       call cebcx(local_scratch,n3d,vec_z,n3d,work_z,maxvec+1,n3d,it+1,1)
       title='Final Value of Exp(-i*H*t) on Initial Vector Final'
       call prntcmn(title,local_scratch,n3d,1,n3d,1,iout,'e')
       IF(log_iterative(10).or.log_iterative(9)) THEN     
          title='Final Solution'
          call prntcmn(title,local_scratch,n3d,1,n3d,1,iout,'e')
       END IF 
       RETURN
  END IF
  IF (it == 0) THEN
       lanczos_tri(it) = work_z(it)
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          title='lanczos_tri iteration = '//itoc(it)
          call prntcmn(title,lanczos_tri,it+1,1,maxvec+1,1,iout,'e')
       END IF
  ELSE
      IF(log_iterative(10).or.log_iterative(8)) THEN     
         title='lanczos_tri iteration = '//itoc(it)
         call prntcmn(title,lanczos_tri,it+1,1,maxvec+1,1,iout,'e')
         title='work iteration = '//itoc(it)
         call prntcmn(title,work_z,it+1,1,maxvec+1,1,iout,'e')
      END IF
      lanczos_tri(0:it) = lanczos_tri(0:it) - work_z(0:it)
      IF(log_iterative(10).or.log_iterative(8)) THEN     
         title='lanczos_tri iteration = '//itoc(it)
         call prntcmn(title,lanczos_tri,it+1,1,maxvec+1,1,iout,'e')
      END IF
  END IF
  wfn_tst = sqrt ( cdotc( it+1 , lanczos_tri , 1 ,lanczos_tri , 1 ) )
  write(iout,1) wfn_tst
  lanczos_tri(0:it) = work_z(0:it)
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
     call cebcx(local_scratch,n3d,vec_z,n3d,work_z,maxvec+1,n3d,it+1,1)
     IF(log_iterative(10).or.log_iterative(9)) THEN     
        title='Final Solution'
        call prntcmn(title,local_scratch,n3d,1,n3d,1,iout,'e')
     END IF
     RETURN
  END IF
1 FORMAT(/,10x,'RMS = sqrt(< ( v_{i} - v_{i-1}) |( v_{i} - v_{i-1}) >) = ',e15.8)
END SUBROUTINE preconditioned_lanczos_convergence_test_z
!***********************************************************************
!***********************************************************************
           END MODULE Lanczos_Module
!***********************************************************************
!***********************************************************************


