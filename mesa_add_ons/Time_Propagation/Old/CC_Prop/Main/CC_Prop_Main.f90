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
  USE dvrprop_global
  USE Propagation_Module
  USE Hamiltonian_Module
  USE Channel_Module
  USE CC_Prop_Regional_Module
  USE CC_Prop_Module
!
  IMPLICIT NONE
  REAL*4                                   :: secnds
  REAL*4, DIMENSION(10)                    :: del_t
  INTEGER                                  :: iostat 
  CHARACTER (LEN=8)                        :: hamiltonian_source
  CHARACTER (LEN=80)                       :: preconditioner
  CHARACTER (LEN=80)                       :: ham_type
  CHARACTER (LEN=80)                       :: type_storage
!
! Read in the Input
!
  time(1)=secnds(0.0)
  CALL CC_Prop_Input(hamiltonian_source,ham_type,preconditioner,type_storage)
  time(2)=secnds(0.0)
  del_t(1) = time(2) - time(1)
  WRITE(iout,1)
  WRITE(iout,2) del_t(1)
  WRITE(iout,1)
!
!        Branch depending on whether the Hamiltonian is generated using the FEDVR library
!        or from another source.
!
                                                                       Hamiltonian_Construction_Section:              &
!
  IF(hamiltonian_source == 'internal') THEN
!
!    Calculate the Global FEDVR Matrices
!
     CALL CC_Prop_Space
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
!             and for the Potential.  Both the Lanczos and Split Operator need these
!             arrays.
!
  ELSE IF(hamiltonian_source == 'external') THEN 
     CALL Read_Ham(ham_type,type_storage)
     time(3)=secnds(0.0)
     del_t(2) = time(3) - time(2)
     WRITE(iout,1)
     WRITE(iout,4) del_t(2)
     WRITE(iout,1)
     IF(diagonalize_only) THEN
        stop
     END IF
  END IF                                                               Hamiltonian_Construction_Section  
!
!
                                                                       Preconditioner_Construction_section:        &
!
  IF (preconditioner == 'cholesky') THEN
!
                                                                       Hamiltonian_Type_I:         &
!
    IF(ham_type == 'real') THEN
       ALLOCATE(upper_d(n3d*(n3d+1)/2), lower_d(n3d*(n3d+1)/2))
!
                                                                       Storage_Type_I:             &
!
       IF( type_storage == 'triangle' ) THEN
           upper_d = triangle_overlap_d
       ELSE IF( type_storage == 'full' ) THEN
           Call Full_to_Upper(overlap_d,upper_d)
!
       END IF                                                          Storage_Type_I 
!
                                                                       Packed_Test_I:              &
!
       IF(packed) THEN
          ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                               &
                   mat_var(1)%row_index(lenbuf),                                                   &
                   mat_var(1)%packed_columns_d(lenbuf) ,                                           &
                   mat_var(1)%matrix_diagonal_d(n3d) )
          ALLOCATE(mat_var(2)%non_zero_columns (n3d),                                              &
                   mat_var(2)%row_index(lenbuf),                                                   &
                   mat_var(2)%packed_columns_d(lenbuf))
!
       END IF                                                          Packed_Test_I
!
       CALL Cholesky(upper_d,lower_d,print_parameter)
                                                                       Packed_Test_II:             &
       IF(packed) THEN                         
          DEALLOCATE(upper_d, lower_d)
          ALLOCATE(mat_var(4)%non_zero_columns(n3d),                                               &
                   mat_var(4)%row_index(lenbuf),                                                   &
                   mat_var(4)%packed_columns_d(lenbuf),                                            &
                   mat_var(4)%matrix_diagonal_d(n3d) )
!
                                                                       Storage_Type_II:            &
          IF (type_storage == 'triangle') THEN
               Call h_pack_upper(triangle_overlap_d,                                               &
                                 mat_var(4)%matrix_diagonal_d,                                     &
                                 mat_var(4)%packed_columns_d,                                      &
                                 mat_var(4)%non_zero_columns,                                      &
                                 mat_var(4)%row_index,                                             &
                                 mat_var(4)%number,n3d,'overlap',                                  &
                                 print_packed)
               DEALLOCATE(triangle_overlap_d)
          ELSE IF(type_storage == 'full') THEN
               Call pack_matrix(overlap_d,                                                         &
                                mat_var(4)%matrix_diagonal_d,                                      &
                                mat_var(4)%packed_columns_d,                                       &
                                mat_var(4)%non_zero_columns,                                       &
                                mat_var(4)%row_index,                                              &
                                mat_var(4)%number,n3d,'overlap',                                   &
                                print_packed)
               DEALLOCATE(overlap_d)
          END IF                                                       Storage_Type_II     
!
       END IF                                                          Packed_Test_II
!
    ELSE IF(ham_type == 'complex') THEN
       ALLOCATE(upper_z(n3d*(n3d+1)/2), lower_z(n3d*(n3d+1)/2))
!
! 
                                                                       Storage_Type_III:           &
!
       IF( type_storage == 'triangle' ) THEN
           upper_z = triangle_overlap_z
       ELSE IF( type_storage == 'full' ) THEN
           Call Full_to_Upper(overlap_z,upper_z)
       END IF                                                          Storage_Type_III
!
                                                                       Packed_Test_III:            &
       IF (packed) THEN
           ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                              &
                    mat_var(1)%row_index(lenbuf),                                                  &
                    mat_var(1)%packed_columns_z(lenbuf),                                           &
                    mat_var(1)%matrix_diagonal_z(n3d) )
          ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                               &
                   mat_var(2)%row_index(lenbuf),                                                   &
                   mat_var(2)%packed_columns_z(lenbuf))
       END IF                                                          Packed_Test_III 
!
       CALL Cholesky(upper_z,lower_z,print_parameter)
                                                                       Packed_Test_IV:             &
       IF(packed) THEN                         
          DEALLOCATE(upper_z, lower_z)
          ALLOCATE(mat_var(4)%non_zero_columns(n3d),                                               &
                   mat_var(4)%row_index(lenbuf),                                                   &
                   mat_var(4)%packed_columns_z(lenbuf),                                            &
                   mat_var(4)%matrix_diagonal_z(n3d) )
                                                                       Storage_Type_IV:            &
!
          IF (type_storage == 'triangle') THEN
               Call h_pack_upper(triangle_overlap_z,                                               &
                                 mat_var(4)%matrix_diagonal_z,                                     &
                                 mat_var(4)%packed_columns_z,                                      &
                                 mat_var(4)%non_zero_columns,                                      &
                                 mat_var(4)%row_index,                                             &
                                 mat_var(4)%number,n3d,'overlap',                                  &
                                 print_packed)
               DEALLOCATE(triangle_overlap_z)
          ELSE IF(type_storage == 'full') THEN
               Call pack_matrix(overlap_z,                                                         &
                                mat_var(4)%matrix_diagonal_z,                                      &
                                mat_var(4)%packed_columns_z,                                       &
                                mat_var(4)%non_zero_columns,                                       &
                                mat_var(4)%row_index,                                              &
                                mat_var(4)%number,n3d,'overlap',                                   &
                                print_packed)
               DEALLOCATE(overlap_z)
          END IF                                                       Storage_Type_IV     
!
       END IF                                                          Packed_Test_IV
!

    END IF                                                             Hamiltonian_Type_I
!
  END IF                                                               Preconditioner_Construction_section
  
                                                                       Hamiltonian_Pack_Section:  &
  IF (packed) THEN
!
                                                                       Hamiltonian_Type_II:       &
!
!
      IF (ham_type == 'real') THEN
          ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                              &
                   mat_var(3)%row_index(lenbuf),                                                  &
                   mat_var(3)%packed_columns_d(lenbuf),                                           &
                   mat_var(3)%matrix_diagonal_d(n3d) )
!
                                                                       Storage_Type_V:            &
!
          IF( type_storage == 'triangle' ) THEN
              Call h_pack_upper(triangle_hamiltonian_d,                                           &
                                mat_var(3)%matrix_diagonal_d,                                     &
                                mat_var(3)%packed_columns_d,                                      &
                                mat_var(3)%non_zero_columns,                                      &
                                mat_var(3)%row_index,                                             &
                                mat_var(3)%number,n3d,'hamiltonian',                              &
                                print_packed)
              DEALLOCATE(triangle_hamiltonian_d)
          ELSE IF( type_storage == 'full' ) THEN
              Call pack_matrix(hamiltonian_d,                                                     &
                               mat_var(3)%matrix_diagonal_d,                                      &
                               mat_var(3)%packed_columns_d,                                       &
                               mat_var(3)%non_zero_columns,                                       &
                               mat_var(3)%row_index,                                              &
                               mat_var(3)%number,n3d,'hamiltonian',                               &
                               print_packed)
              DEALLOCATE(hamiltonian_d)
!
          END IF                                                       Storage_Type_V
!
      ELSE IF (ham_type == 'complex') THEN
          ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                              &
                   mat_var(3)%row_index(lenbuf),                                                  &
                   mat_var(3)%packed_columns_z(lenbuf),                                           &
                   mat_var(3)%matrix_diagonal_z(n3d) )                  
!
                                                                      Storage_Type_VI:            &
!
          IF( type_storage == 'triangle' ) THEN
              Call h_pack_upper(triangle_hamiltonian_z,                                           &
                                mat_var(3)%matrix_diagonal_z,                                     &
                                mat_var(3)%packed_columns_z,                                      &
                                mat_var(3)%non_zero_columns,                                      &
                                mat_var(3)%row_index,                                             &
                                mat_var(3)%number,n3d,'hamiltonian',                              &
                                print_packed)
              DEALLOCATE(triangle_hamiltonian_z)
          ELSE IF( type_storage == 'full' ) THEN
              Call pack_matrix(hamiltonian_z,                                                     &
                               mat_var(3)%matrix_diagonal_z,                                      &
                               mat_var(3)%packed_columns_z,                                       &
                               mat_var(3)%non_zero_columns,                                       &
                               mat_var(3)%row_index,                                              &
                               mat_var(3)%number,n3d,'hamiltonian',                               &
                               print_packed)
              DEALLOCATE(hamiltonian_z)
!
          END IF                                                       Storage_Type_VI
!
      END IF                                                           Hamiltonian_Type_II                     
!
  END IF                                                               Hamiltonian_Pack_Section
!
!    Here is where we have to begin to parallelize the code.
!    Arrays have to be allocated to different processors and I am not sure of
!    the best way to proceed.
!
  time(1)=secnds(0.0)!
  OPEN(UNIT=99,FILE='initial-wavefunction',ACCESS='sequential',                      &
       FORM='unformatted',IOSTAT=IOSTAT,STATUS='unknown')
  CALL Iterative_Data
!        
!    The iterative methods need lots of memory in order to be efficient.
!    The big arrays are n3d*maxvec in size where maxvec is the maximum
!    method is much less demanding and this is a big plus for many applications.
!      
  CALL Propagation_Driver
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
4 FORMAT(/,10X,'Time to Read, Factor and Diagonalize Matrices = ',f15.8)
5 FORMAT(/,10x,'Total time for the Iterative Propagation       = ',f15.8)
  stop
END PROGRAM CC_Prop_Main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
