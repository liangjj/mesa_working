!***********************************************************************
! Lanczos_Module
!**begin prologue     Lanczos_Module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Lanczos_Library
!**purpose            Contains all of the major subroutines 
!***                  to either find the eigenvalues of a 
!***                  real symmetric matrix or to solve a set of linear
!***                  using the Lanczos algorithm.
!***                  Explicit interfaces are used to allow
!***                  a transparent use of generic subroutines which work
!***                  for both real and complex vectors.  
!***description       Given a starting vector, a number of iterations
!***                  are performed until the eigenvalues or solution
!***                  converge to a fixed accuracy criterion.
!**references
!**modules needed     See USE statements below
!**comments           In this portable version I have disabled all unnecessary
!**                   writing to files.  The original Fortran is commented out.
!**end prologue       Lanczos_Module
!***********************************************************************
!***********************************************************************
                           MODULE Lanczos_Module
                           USE full_matrix_vector_multiply_module
                           USE packed_matrix_vector_multiply_module
                           USE Iterative_Global
                           USE Pack_Global
                           USE input_output
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
                           INTERFACE lanczos_solve
             MODULE PROCEDURE lanczos_solve_d,                     &
                              lanczos_solve_z               
                       END INTERFACE lanczos_solve
!
                           INTERFACE lanczos_convergence_test
             MODULE PROCEDURE lanczos_convergence_test_d,          &
                              lanczos_convergence_test_z  
                       END INTERFACE lanczos_convergence_test
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
  LOGICAL                       :: dollar
  LOGICAL                       :: logkey
  LOGICAL                       :: linsys
  CHARACTER (LEN=80)            :: chrkey
  REAL*8                        :: fpkey
  INTEGER                       :: intkey 
!
! Read in the data
!
  IF ( dollar('$iterative_data',card,cpass,inp) ) then
      iterative_method=chrkey(card,'iterative_method','lanczos',' ')
      linsys=logkey(card,'linear_system',.false.,' ')
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
!**end prologue       lanczos_d
!***********************************************************************
  SUBROUTINE lanczos_d(wave_function)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                :: wave_function
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
  anorm = 1.d0 / sqrt (ddot(n3d , wave_function(:) , 1, wave_function(:) , 1 ))
  vec_d(:,0) = anorm * wave_function(:)
  DO i=0, maxvec
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
        call prntfmn(title,vec_d(:,i),n3d,1,n3d,maxvec+1,iout,'e')
     END IF
!
     call lanczos_vec(wave_function,vec_d,i)
!
!
!           If the number of iterations has reached the size of the
!           original matrix without breakdown, we have to be converged.
!
    IF( convergence == 'maximal') THEN
        return
!
!           We are actually converged to the desired tolerance.
!
    ELSE IF(convergence == 'converged') THEN
        WRITE(iout,5) i  
        return
    END IF    
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
  SUBROUTINE lanczos_z(wave_function)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)            :: wave_function
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
  anorm = 1.d0 / sqrt (cdotc(n3d , wave_function(:) , 1, wave_function(:) , 1 ))
  vec_z(:,0) = anorm * wave_function(:)
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
    call lanczos_vec(wave_function,vec_z,i)
    IF( convergence == 'maximal') THEN
        WRITE(iout,4)
         return
    ELSE IF(convergence == 'converged') THEN
        WRITE(iout,5) i  
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
  SUBROUTINE lanczos_vec_d(wave_function,lanczos_basis,it)
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d,0:maxvec)                :: lanczos_basis
  REAL*8, DIMENSION(n3d)                         :: wave_function
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
  call column_packed_matrix_on_vector (lanczos_basis(1:n3d,it), h_vec_d,      &
                                       matrix_diagonal_d,                     &
                                       packed_columns_v_d,                    &
                                       non_zero_columns_v,row_index_v,n3d)
  IF(log_iterative(10).or.log_iterative(3)) THEN
     title='H_on_Lanczos_Vector iteration = '//itoc(it)
     call prntfmn(title,h_vec_d,n3d,1,n3d,1,iout,'e')
  END IF
!
!
!             Solve Linear System of the projected Problem
!
  call lanczos_solve(h_mat_d,it)
!
!             Test convergence
!
! 
  call lanczos_convergence_test(wave_function,it)
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
              title='before any reorthogonalization'
              call prntfmn(title,lanczos_basis,n3d,it+2,n3d,maxvec+1,iout,'e')
          END IF
!
!           Performing a reorthogonalize orthogonalization can help
!           prevent linear dependence.  Its an option if needed.
!
          IF(orthogonalize /= 'none') THEN        
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
             IF(log_iterative(10)) THEN
                title='after reorthogonalization'
                call prntfmn(title,lanczos_basis,n3d,it+2,n3d,maxvec+1,iout,'e')
             END IF
          END IF
      END IF
!
!          Compute b(it+1)  b(i) = <V_i|V_i> 
!
      b(it+1) = SQRT( ddot(n3d , lanczos_basis(:,it+1) , 1 ,       &
                                 lanczos_basis(:,it+1) , 1) ) 
      IF(log_iterative(10).or.log_iterative(2)) THEN
         title='Lanczos Beta iteration = '//itoc(it)
         call prntfmn(title,b(it+1),1,1,maxvec+1,maxvec+1,iout,'e')
      END IF
!
!          Test for breakdown.  It occurs if b(i) becomes zero.
!
      IF( b(it+1) <= eps ) THEN
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
            call prntfmn(title,b(it+1),1,1,maxvec+1,maxvec+1,iout,'e')
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
  SUBROUTINE lanczos_vec_z(wave_function,lanczos_basis,it)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d,0:maxvec)            :: lanczos_basis
  COMPLEX*16, DIMENSION(n3d)                     :: wave_function
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
  call column_packed_matrix_on_vector (lanczos_basis(1:n3d,it), h_vec_z,      &
                                       matrix_diagonal_z, packed_columns_v_z, &
                                       non_zero_columns_v,row_index_v,n3d)
  IF(log_iterative(10).or.log_iterative(3)) THEN
     title='H_on_Lanczos_Vector iteration = '//itoc(it)
     call prntcmn(title,h_vec_z,n3d,1,n3d,1,iout,'e')
  END IF
!
!             Solve Linear System of the projected Problem.
!
  call lanczos_solve(h_mat_z,it)
!
!            Test convergence
!
  call lanczos_convergence_test(wave_function,it)
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
             title='before any reorthogonalization'
             call prntcmn(title,lanczos_basis,n3d,it+2,n3d,maxvec+1,iout,'e')
          END IF
!
!           Performing a reorthogonalize orthogonalization can help
!           prevent linear dependence.  Its an option if needed.
!
          IF(orthogonalize /= 'none') THEN        
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
                call prntfmn(title,lanczos_basis,n3d,it+2,n3d,maxvec+1,iout,'e')
             END IF
          END IF
      END IF
!
!          Compute b(it+1)
!
      b(it+1) = SQRT( cdotc(n3d , lanczos_basis(:,it+1) , 1 ,       &
                              
                                  lanczos_basis(:,it+1) , 1) ) 
      IF(log_iterative(10).or.log_iterative(2)) THEN
         title='Lanczos Beta iteration = '//itoc(it)
         call prntfmn(title,b(it+1),1,1,maxvec+1,maxvec+1,iout,'e')
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
!**begin prologue     lanczos_solve_d
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
!**end prologue       lanczos_solve_d
!***********************************************************************
  SUBROUTINE lanczos_solve_d(h,it)
  IMPLICIT NONE
  REAL*8, DIMENSION(0:maxvec,0:2)                :: h
  INTEGER                                        :: i, it
  REAL*8                                         :: ddot
  CHARACTER(LEN=4)                               :: itoc
!
!
!  Now that we have this we can compute the value of a(it) as the
!
!                  a(i) = <V_i|h|V_i> 
!
  a(it) =  ddot(n3d , vec_d(:,it) , 1 , h_vec_d(:) , 1 )
  IF(log_iterative(10).or.log_iterative(1)) THEN
     title='Lanczos Alpha iteration = '//itoc(it)
     call prntfm(title,a(it),1,1,maxvec+1,1,iout)
  END IF
!
!
  rhs_tri_d(it) = ddot(n3d , vec_d(:,it) , 1 , rhs_d(:) , 1 )
  soln_tri_d(it) = rhs_tri_d(it)
!
! Zero Hamiltonian
!
  h(0:it,0:2) = 0.d0
!
! Fill the Hamiltonian
!
  DO i=0,it
     h(i,0) = a(i)
  END DO
  DO i=0,it-1
     h(i,1) = b(i+1)
     soln_tri_d(i) = rhs_tri_d(i)
  END DO
  h(0:it-1,2) = h(0:it-1,1)    
  IF(log_iterative(10).or.log_iterative(5)) THEN
     title = 'Lanczos Diagonal iteration = '//itoc(it)
     call prntfm(title,h(:,0),it+1,1,maxvec+1,1,iout)
     title = 'Lanczos Off-Diagonal iteration = '//itoc(it)
     call prntfm(title,h(:,1),it+1,1,maxvec+1,1,iout)
  END IF
!
!     Lets Solve the Equations.
!
  call dgtsv(it+1,1,h(:,1),h(:,0),h(:,2),soln_tri_d,it+1,info)
  title='Solution in Lanczos Basis iteration = '//itoc(it)
  call prntfm(title,soln_tri_d,it+1,1,maxvec+1,1,iout)
!
END SUBROUTINE lanczos_solve_d
!***********************************************************************
!***********************************************************************
!**begin prologue     lanczos_solve_z
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
!**end prologue       lanczos_solve_z
!***********************************************************************
  SUBROUTINE lanczos_solve_z(h,it)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(0:maxvec,0:2)                :: h
  INTEGER                                            :: i, it
  COMPLEX*16                                         :: cdotc
  CHARACTER(LEN=4)                                   :: itoc
!
!            Now that we have this we can compute the value of a(it) as the
!                                 a(i) = <V_i|h|V_i> 
!
!
  a(it) =  cdotc(n3d , vec_z(:,it) , 1 , h_vec_z(:) , 1 )
  IF(log_iterative(10).or.log_iterative(1)) THEN
     title='Lanczos Alpha iteration = '//itoc(it)
     call prntfm(title,a(it),1,1,maxvec+1,1,iout)
  END IF
!
!
  rhs_tri_z(it) = cdotc(n3d , vec_z(:,it) , 1 , rhs_z(:) , 1 )   
  soln_tri_z(it) = rhs_tri_z(it)
!
!
! Zero Hamiltonian
!
  h(0:it,0:2) = 0.d0
!
! Fill the Hamiltonian
!
  DO i=0,it
     h(i,0) = a(i)
  END DO
  DO i=0,it-1
     h(i,1) = b(i+1)
     soln_tri_z(it) = rhs_tri_z(it)
  END DO
  h(0:it-1,2) = h(0:it-1,1)    
  IF(log_iterative(10).or.log_iterative(5)) THEN
     title = 'Lanczos Diagonal iteration = '//itoc(it)
     call prntcm(title,h(:,0),it+1,1,maxvec+1,1,iout)
     title = 'Lanczos Off-Diagonal iteration = '//itoc(it)
     call prntcm(title,h(:,1),it+1,1,maxvec+1,1,iout)
  END IF
!
!    Lets Solve the Equations.
!
  call zgtsv(it+1,1,h(:,1),h(:,0),h(:,2),soln_tri_z,it+1,info)
  title='Solution in Lanczos Basis iteration = '//itoc(it)  
  call prntcm(title,soln_tri_z,it+1,1,maxvec+1,1,iout)
!
END SUBROUTINE lanczos_solve_z
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
  SUBROUTINE lanczos_convergence_test_d(wave_function,it)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                           :: wave_function
  INTEGER                                        :: i, it
  REAL*8                                         :: ddot
  CHARACTER(LEN=4)                               :: itoc
!
!               
!
  IF (it == 0) THEN
       lanczos_tri_d(it) = soln_tri_d(it)
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          title='lanczos_tri iteration = '//itoc(it)
          call prntfm(title,lanczos_tri_d,it+1,1,maxvec+1,1,iout)
       END IF
       WRITE(iout,1)
       RETURN
  ELSE IF( it+1 == n3d) THEN
!
!             We have maximized the number of iterations, we are totally finished.  
       convergence = 'maximal'
       call ebcx(wave_function,n3d,vec_d,n3d,soln_tri_d,maxvec+1,n3d,it+1,1)
       title='Final Solution'
       call prntfm(title,wave_function,n3d,1,n3d,1,iout)
       WRITE(iout,2)
       RETURN 
  ELSE
       work_d(0:it) = soln_tri_d(0:it) - lanczos_tri_d(0:it-1)
       wfn_tst = sqrt ( ddot( it+1 , work_d , 1 , work_d , 1 ) )
       write(iout,3) wfn_tst
       lanczos_tri_d(0:it) = soln_tri_d(0:it)
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          title='work iteration = '//itoc(it)
          call prntfm(title,work_d,it+1,1,maxvec+1,1,iout)
       END IF
!
!             We are below the maximal number of iterations.  Test convergence.
!
       IF(wfn_tst <= cnverg) THEN
!
!                  
!
          convergence = 'converged'
          call ebcx(wave_function,n3d,vec_d,n3d,soln_tri_d,maxvec+1,n3d,it+1,1)
          title='Final Solution'
          call prntfm(title,wave_function,n3d,1,n3d,1,iout)
       END IF
       RETURN
  END IF
1 FORMAT(/,10x,'Iteration = 0 No Test')
2 FORMAT(/,10x,'Iteration = n No Test')
3 FORMAT(/,10x,'RMS = sqrt(< ( v_{i} - v_{i-1}) |( v_{i} - v_{i-1}) >) = ',e15.8)
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
!******************************************************************
  SUBROUTINE lanczos_convergence_test_z(wave_function,it)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)                           :: wave_function
  INTEGER                                            :: i, it
  COMPLEX*16                                         :: cdotc
  CHARACTER(LEN=4)                                   :: itoc
  IF (it == 0) THEN
       lanczos_tri_z(it) = soln_tri_z(it)
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          title='lanczos_tri iteration = '//itoc(it)
          call prntcm(title,lanczos_tri_z,it+1,1,maxvec+1,1,iout)
       END IF
       WRITE(iout,1)
       RETURN
  ELSE IF( it+1 == n3d) THEN
!
!             We have maximized the number of iterations, we are totally finished.  
       convergence = 'maximal'
       call cebcx(wave_function,n3d,vec_d,n3d,soln_tri_d,maxvec+1,n3d,it+1,1)
       title='Final Solution'
       call prntcm(title,wave_function,n3d,1,n3d,1,iout)
       WRITE(iout,2) 
       RETURN
  ELSE
       work_z(0:it) = soln_tri_z(0:it) - lanczos_tri_z(0:it-1)
       IF(log_iterative(10).or.log_iterative(8)) THEN     
          title='lanczos_tri iteration = '//itoc(it)
          call prntfm(title,lanczos_tri_z,it+1,1,maxvec+1,1,iout)
          title='work iteration = '//itoc(it)
          call prntfm(title,work_z,it+1,1,maxvec+1,1,iout)
       END IF
       wfn_tst = sqrt ( cdotc( it+1 , work_z , 1 , work_z , 1 ) )
       write(iout,3) wfn_tst
       lanczos_tri_z(0:it) = soln_tri_z(0:it)
!
!             We are below the maximal number of iterations.  Test convergence.
!
       IF(wfn_tst <= cnverg) THEN
!
!                  
!
          convergence = 'converged'
          call cebcx(wave_function,n3d,vec_z,n3d,soln_tri_z,maxvec+1,n3d,it+1,1)
          title='Final Solution'
          call prntcm(title,wave_function,n3d,1,n3d,1,iout)
        END IF
        RETURN
   END IF
1 FORMAT(/,10x,'Iteration = 0 No Test')
2 FORMAT(/,10x,'Iteration = n No Test')
3 FORMAT(/,10x,'RMS = sqrt(< ( v_{i} - v_{i-1}) |( v_{i} - v_{i-1}) >) = ',e15.8)
!
END SUBROUTINE lanczos_convergence_test_z
!***********************************************************************
!***********************************************************************
END MODULE Lanczos_Module
!***********************************************************************
!***********************************************************************



