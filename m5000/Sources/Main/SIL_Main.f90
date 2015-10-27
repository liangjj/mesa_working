!deck SIL_Main
!**begin prologue     SIL_Main
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            main propagator program
!**description        calls routines to input data, construct FEDVR propagators
!**                   for a close-coupled wave function.
!**                   
!**references
!**routines called      Name                    Location    
!                       ----                    --------
!                     Input_Prop             SIL_Sources
!                     Space_Prop             SIL_Sources
!                     Regional_Matrices      SIL_Modules:regional_module
!                     v_couple               SIL_Sources
!                     Propagation_Driver     SIL_Modules:Propagation_Module

!
!**modules used         Name                    Location               
!                       ----                    --------
!                    dvrprop_global          Modules   
!                    Propagation_Module      SIL_Modules
!
!**end prologue       SIL_Main
  PROGRAM SIL_Main
  USE Propagation_Module
  USE Data_Subroutines_Module
  USE Packed_Matrix_Module
  USE Atomic_State_Matrix_Module
!
  IMPLICIT NONE
  REAL*4                                   :: secnds
  REAL*4, DIMENSION(10)                    :: del_t
  INTEGER                                  :: iostat 
  INTEGER                                  :: i 
  INTEGER                                  :: lenh 
  INTEGER                                  :: lenth
  CHARACTER (LEN=80)                       :: preconditioner
  CHARACTER (LEN=80)                       :: ham_type
  CHARACTER (LEN=80)                       :: home_directory
!
!
! Read in the Input
!
  Call GetEnv('SOURCES',home_directory)
  lenh=lenth(home_directory)
  home_directory=home_directory(1:lenh)//'/SIL_Propagator'
  nunits=10
  n_dir=6
  ALLOCATE(names(1:nunits),len_name(1:nunits),File_Directory(1:n_dir),len_dir(1:n_dir))
  Call Command_Line(home_directory=home_directory)
  time(1)=secnds(0.0)
  CALL SIL_Input(ham_type,preconditioner)                 ! This routine is located in the Data_Subroutines_Module
                                                          ! It sets up a broad array of variables needed for the
                                                          ! including where to find certain quantities that need 
                                                          ! to read from disk and the manner in which various 
                                                          ! arrays are formatted.  The user should look at that 
                                                          ! routine in detail to see what is going on.
  time(2)=secnds(0.0)
  del_t(1) = time(2) - time(1)
  WRITE(iout,1)
  WRITE(iout,2) del_t(1)
  WRITE(iout,1)
!
!        Branch depending on whether the Hamiltonian is generated using the 
!        FEDVR library or from another source.
!
  WRITE(iout,*) 'Opening File eigenstates to hold data'  ! Open a file to store data from the calculation
                                                         ! If one needs to restart the calculation the wavefunction
                                                         ! from the previous timestep needs to be dumped.
!
  OPEN(UNIT=20,FILE=File_Directory(4)(1:len_dir(4))//'/eigenstates.dat',ACCESS='sequential',   &
       FORM='unformatted',IOSTAT=IOSTAT,STATUS='unknown')
!
!
!
!    ALLOCATE the Memory for the Matrices used to Store the Wavefunction
!    and for the Potential.
!
  CALL Directives(ham_type)                ! This routine reads options and data required in the calculation.
                                           ! The data includes matrix directives, printing, buffers and formatting
                                           ! control.
  time(3)=secnds(0.0)
  del_t(2) = time(3) - time(2)
  WRITE(iout,1)
  WRITE(iout,3) del_t(2)
!
!
  time(1) = secnds(0.0)
  IF(input_matrices(1:6) == 'packed'.or.output_matrices(1:6) == 'packed')  THEN
!
     Call Packed_Matrix(ham_type)                ! The input/output matrices are/will be stored in a packed format
                                                 ! with no zeros but also with no recognition that there is angular
                                                 ! symmetry.  The routines used are in the packed_matrix module.
  ELSE IF(input_matrices(1:22) == 'using_angular_symmetry') THEN
     Call Pack_State_Matrices                    ! The matrices are also packed as above but using symmetry.
                                                 ! The routines used are in the atomic_state_matrix module.
  END IF
  IF(diagonalize_only) THEN
     stop
  END IF
!
!    Here is where we have to begin to parallelize the code.
!    Arrays have to be allocated to different processors and I am not sure of
!    the best way to proceed.
!
  time(1)=secnds(0.0)
  OPEN(UNIT=99,FILE=File_Directory(4)(1:len_dir(4))//'/initial-wavefunction.dat',          &
       ACCESS='sequential',FORM='unformatted',IOSTAT=IOSTAT,STATUS='unknown')
  CALL Iterative_Data                           ! This routine inputs the data required for the Lanczos iteration
                                                ! such as number of iterations, allowed vector space size, 
                                                ! convergence, etc.
!       
!    The iterative methods need lots of memory in order to be efficient.
!    The big arrays are n3d*maxvec in size where maxvec is the maximum size of the vector space.
!
!                     On to the actual propagation.
  CALL Propagation_Driver                       ! This is the propagation driver which does most of the 
                                                ! allocation of storage and then calls the routines which
                                                ! do the propagation.  The code caan deal with the ordinary
                                                ! Lanczos with unit metric as well as a generalized problem where
                                                ! one uses the Cholesky decomposition of the overlap either to
                                                ! transform the Hamiltonian explicity so that it is in standard
                                                ! form or directly, using the Cholesky decomposition at each
                                                ! to solve for the next iteration vector.  In order to use the
                                                ! the direct method one needs to solve two sets of triangular
                                                ! linear equations at each iteration.  Since the matrix has been
                                                ! factored, the number of operations to perform the solution is
                                                ! the same as the matrix-vector multiplication by a symmetric
                                                ! matrix.  This can be much less work than the explicit approach
                                                ! but a lot depends on the size of the problem, symmetry, zeros and
                                                ! a lot of other considerations.  It needs to be looked at carefully
                                                ! in each specific case.
  CLOSE(20)
  time(2) = secnds(0.0)
  del_t(1) = time(2) - time(1)
  WRITE(iout,1)
  WRITE(iout,4) del_t(1)
  WRITE(iout,1)
1 FORMAT('***********************************************'                           &
         '*************************')
2 FORMAT(/,10X,'Time to Input Basic Data               = ',f15.8)
3 FORMAT(/,10X,'Time to Read in Data = ',f15.8)
4 FORMAT(/,10x,'Time for the Iterative Propagation       = ',f15.8)
  stop
END PROGRAM SIL_Main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
