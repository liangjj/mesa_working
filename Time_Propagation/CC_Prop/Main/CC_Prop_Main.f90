!deck CC_Prop_Main
!**begin prologue     CC_Prop_Main
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
!                     Input_Prop             CC_Prop_Sources
!                     Space_Prop             CC_Prop_Sources
!                     Regional_Matrices      CC_Prop_Modules:regional_module
!                     v_couple               CC_Prop_Sources
!                     Propagation_Driver     CC_Prop_Modules:Propagation_Module

!
!**modules used         Name                    Location               
!                       ----                    --------
!                    dvrprop_global          Modules   
!                    Propagation_Module      CC_Prop_Modules
!
!**end prologue       CC_Prop_Main
  PROGRAM CC_Prop_Main
  USE Propagation_Module
  USE CC_Prop_Regional_Module
  USE Data_Subroutines_Module
  USE Non_Packed_Matrix_Module
  USE Packed_Matrix_Module
  USE Atomic_State_Matrix_Module
!
  IMPLICIT NONE
  REAL*4                                   :: secnds
  REAL*4, DIMENSION(10)                    :: del_t
  INTEGER                                  :: iostat 
  INTEGER                                  :: lenth
  INTEGER                                  :: len 
  CHARACTER (LEN=80)                       :: hamiltonian_source
  CHARACTER (LEN=80)                       :: preconditioner
  CHARACTER (LEN=80)                       :: ham_type
  CHARACTER (LEN=80)                       :: type_storage
  CHARACTER (LEN=80)                       :: file_directory

!
! Read in the Input
!

  time(1)=secnds(0.0)
  CALL CC_Prop_Input(hamiltonian_source,ham_type,preconditioner,type_storage,file_directory)
  len=lenth(file_directory)
  time(2)=secnds(0.0)
  del_t(1) = time(2) - time(1)
  WRITE(iout,1)
  WRITE(iout,2) del_t(1)
  WRITE(iout,1)
!
!        Branch depending on whether the Hamiltonian is generated using the 
!        FEDVR library or from another source.
!
  WRITE(iout,*) 'Opening File eigenstates to hold data'
  OPEN(UNIT=20,FILE=file_directory(1:len)//'/eigenstates.dat',ACCESS='sequential', FORM='unformatted',   &
       IOSTAT=IOSTAT,STATUS='unknown')
!
                      Hamiltonian_Construction_Section:                          &
!
  IF(hamiltonian_source == 'internal') THEN
!
!    Set the FEDVR variables to default values
!
     Call Set_FEDVR_Defaults
! 
!    Calculate the Global FEDVR Matrices
!
     CALL CC_Prop_Space
!
     time(3)=secnds(0.0)
     del_t(2) = time(3) - time(2)
!
!    Compute the Regional Matrices from the FEDVR Global Matrices
!
     CALL CC_Prop_Regional_Matrices
     time(4)=secnds(0.0)
     del_t(3) = time(4)-time(3)
     WRITE(iout,1)
     WRITE(iout,3) del_t(2:3)
     WRITE(iout,1)
!
!    ALLOCATE the Memory for the Matrices used to Store the Wavefunction
!    and for the Potential.  Both the Lanczos and Split Operator need these
!    arrays.
!
  ELSE IF(hamiltonian_source == 'external') THEN 
     CALL Directives(ham_type,file_directory)
     time(3)=secnds(0.0)
     del_t(2) = time(3) - time(2)
     WRITE(iout,1)
     WRITE(iout,4) del_t(2)

!
  END IF              Hamiltonian_Construction_Section  
!
  time(1) = secnds(0.0)
  IF (input_matrices == 'full_matrices') THEN
      Call Small_Matrix(ham_type,file_directory,ham_file)  
      WRITE(iout,1)
      IF(diagonalize_only) THEN
         stop
      END IF
      IF (test_s_inverse) THEN
          stop
      END IF
      WRITE(iout,*) 'Closing File eigenstates to hold data'
      CLOSE(20)
!
  ELSE IF(input_matrices(1:6) == 'packed')  THEN
!
      Call Packed_Matrix(ham_type,file_directory)  
!
  ELSE IF(input_matrices(1:22) == 'using_angular_symmetry') THEN
      Call Pack_State_Matrices(file_directory)
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
  OPEN(UNIT=99,FILE=file_directory(1:len)//'/initial-wavefunction.dat',          &
       ACCESS='sequential',FORM='unformatted',IOSTAT=IOSTAT,STATUS='unknown')
  CALL Iterative_Data(hamiltonian_source)
!        
!    The iterative methods need lots of memory in order to be efficient.
!    The big arrays are n3d*maxvec in size where maxvec is the maximum
!    method is much less demanding and this is a big plus for many applications.
!
!                     On to the actual propagation.
  CALL Propagation_Driver(file_directory)
  CLOSE(20)
  time(2) = secnds(0.0)
  del_t(1) = time(2) - time(1)
  WRITE(iout,1)
  WRITE(iout,5) del_t(1)
  WRITE(iout,1)
1 FORMAT('***********************************************'                           &
         '*************************')
2 FORMAT(/,10X,'Time to Input Basic Data               = ',f15.8)
3 FORMAT(/,10X,'Time to Compute Global FEDVR Spatial ',                              &
         /,10x,'Matrices                               = ',f15.8,/,                  &
         /,10x,'Time to Compute the Regional Spatial '                               &
         /,10x,'Matrices                               = ',f15.8)
4 FORMAT(/,10X,'Time to Read in Data = ',f15.8)
5 FORMAT(/,10x,'Time for the Iterative Propagation       = ',f15.8)
  stop
END PROGRAM CC_Prop_Main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
