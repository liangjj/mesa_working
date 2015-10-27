!***********************************************************************
                           MODULE Cholesky_Decomposition
                           USE dvrprop_global
                           USE CC_Prop_Module
                           USE Preconditioner_Module
                           USE Iterative_Global
                           USE Pack_Global
                           USE Pack_Hamiltonian_Module
!
                           IMPLICIT NONE


  CHARACTER (LEN=80)                       :: ham_file
  CHARACTER (LEN=80)                       :: ham_type
  INTEGER                                  :: number_of_splines
  INTEGER                                  :: number_of_channels
  INTEGER                                  :: number_of_correlation_terms
  INTEGER                                  :: chan_size
  INTEGER                                  :: number_of_splines_i
  INTEGER                                  :: number_of_splines_j
  INTEGER                                  :: first_spline
  INTEGER                                  :: last_spline
  LOGICAL                                  :: ifeig=.false.
  LOGICAL                                  :: diagonalize_only
  LOGICAL                                  :: read_only
  LOGICAL                                  :: channel_matrices
  LOGICAL, DIMENSION(2)                    :: spline_removal
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



                           Contains
!***********************************************************************
!***********************************************************************
!deck read_ham
!***begin prologue     read_ham
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            read input or B-spline Hamiltonian
!***
!***references
!***routines called
!***end prologue       read_ham
!
  SUBROUTINE read_ham
  IMPLICIT NONE
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  LOGICAL                               :: dollar, logkey
  INTEGER                               :: intkey
  INTEGER                               :: i
  REAL*8                                :: fpkey
  INTEGER                               :: len, lenth
  INTEGER                               :: iostat
!
!
  IF ( dollar('$hamiltonian_parameters',card,cpass,inp) ) then
       i0stat=chrkey(card,'initial_state','from_input',' ')
       ham_file=chrkey(card,'hamiltonian_file_name','from_file',' ')
       ham_type=chrkey(card,'hamiltonian_type','real',' ')
       triangle=logkey(card,'triangular_storage',.false.,' ')
       print_parameter=logkey(card,'print=on',.false.,' ')
       print_packed=logkey(card,'print_packed_matrix=on',.false.,' ')
       diagonalize_only=logkey(card,'diagonalize_only',.false.,' ')
       rows_to_print = intkey(card,'rows_to_print',20,' ')
       eigenvectors_to_print = intkey(card,'eigenvectors_to_print',5,' ')
       preconditioner=chrkey(card,'preconditioner','none',' ')
       read_only = logkey(card,'read_only',.false.,' ')
       channel_matrices = logkey(card,'channel_matrices',.false.,' ')
       IF (channel_matrices) THEN
           number_of_channels = intkey(card,'number_of_channels',1,' ')
           number_of_splines = intkey(card,'number_of_splines',1,' ')
           number_of_correlation_terms = intkey(card,'number_of_correlation_terms',  &
                                                0,' ')
       END IF
       len=lenth(ham_file)
       WRITE(iout,*) 'Hamiltonian File Name = ',ham_file(1:len)
       WRITE(iout,*) 'Non Orthogonal Basis = ',non_orth
  END IF
  IF( ham_file == 'from_input') THEN
      write(iout,*) 'Hamiltonian File from Input'
      n3d = intkey(card,'matrix_size',2,' ')
  ELSE 
      write(iout,*) 'Hamiltonian File from Disk'
      OPEN(UNIT=50,FILE=ham_file(1:len),ACCESS='sequential',                             &
           FORM='unformatted',IOSTAT=IOSTAT,STATUS='old')
      READ(50) n3d
      IF (channel_matrices) THEN
          spline_removal(1) = logkey(card,'remove_first_spline',.false.,' ')
          spline_removal(2) = logkey(card,'remove_last_spline',.false.,' ')
          ALLOCATE(channel_mat(number_of_channels,number_of_channels),                   &
                   channel_buf(number_of_channels,number_of_splines))
      END IF
  END IF
  IF (packed) THEN
      WRITE(iout,1)
  END IF
  IF(i0stat /= 'from_disk') THEN
     write(iout,*) 'Opening File eigenstates to hold data'
     OPEN(UNIT=20,FILE='eigenstates',ACCESS='sequential', FORM='unformatted',            &
          IOSTAT=IOSTAT,STATUS='unknown')
  END IF
  IF(in_core) THEN
     lenbuf = n3d*(n3d-1)/2
  ELSE
     lenbuf = min(lenbuf,n3d*(n3d-1)/2)
  END IF
  rows_to_print = min(rows_to_print,n3d)      
  eigenvectors_to_print = min(eigenvectors_to_print,n3d)
  ALLOCATE(rowlab(n3d),collab(n3d))
  DO i=1,n3d
     rowlab(i) = 'Row '//itoc(i)
     collab(i) = 'Col '//itoc(i)
  END DO
  lwork=10*n3d 
  WRITE(iout,2) n3d
  WRITE(iout,*) 'Hamiltonian Type = ',ham_type(1:len)
  IF (ham_type == 'real') THEN
!
!     
      IF (triangle) THEN
!
!     
!         Use Triangular Storage.
!
          ALLOCATE(triangle_hamiltonian_d(n3d*(n3d+1)/2))
          IF(non_orth) THEN      
             ALLOCATE(triangle_overlap_d(n3d*(n3d+1)/2))
          END IF
          Call Compute_Eigenstates(triangle_hamiltonian_d, triangle_overlap_d) 
          IF(preconditioner == 'cholesky_decomposition') THEN
             ALLOCATE(upper_d(n3d*(n3d+1)/2), lower_d(n3d*(n3d+1)/2)) 
             write(iout,3)
             upper_d = triangle_overlap_d
             CALL Cholesky(upper_d,lower_d,print_parameter)
             IF (packed) THEN
!
!                Use packed storage from now on.
!
                 DEALLOCATE(upper_d, lower_d)
                 ALLOCATE(mat_var(3)%non_zero_columns(n3d),                        &
                          mat_var(3)%row_index(lenbuf),                            &
                          mat_var(3)%packed_columns_d(lenbuf),                     &
                          mat_var(3)%matrix_diagonal_d(n3d) )
                 Call h_pack_upper(triangle_hamiltonian_d,                         &
                                  mat_var(3)%matrix_diagonal_d,                    &
                                  mat_var(3)%packed_columns_d,                     &
                                  mat_var(3)%non_zero_columns,                     &
                                  mat_var(3)%row_index,                            &
                                  mat_var(3)%number,n3d,'hamiltonian',             &
                                  print_packed)
                 DEALLOCATE(triangle_hamiltonian_d)
                 ALLOCATE(mat_var(4)%non_zero_columns(n3d),                        &
                          mat_var(4)%row_index(lenbuf),                            &
                          mat_var(4)%packed_columns_d(lenbuf),                     &
                          mat_var(4)%matrix_diagonal_d(n3d) )
                 Call h_pack_upper(triangle_overlap_d,                             &
                                  mat_var(4)%matrix_diagonal_d,                    &
                                  mat_var(4)%packed_columns_d,                     &
                                  mat_var(4)%non_zero_columns,                     &
                                  mat_var(4)%row_index,                            &
                                  mat_var(4)%number,n3d,'overlap',                 &
                                  print_packed)
                 DEALLOCATE(triangle_overlap_d)
             END IF
          ELSE
             IF (packed) THEN
!
!                Use packed storage from now on.
!
                 ALLOCATE(mat_var(3)%non_zero_columns(n3d),                        &
                          mat_var(3)%row_index(lenbuf),                            &
                          mat_var(3)%packed_columns_d(lenbuf),                     &
                          mat_var(3)%matrix_diagonal_d(n3d) )
                 Call h_pack_upper(triangle_hamiltonian_d,                         &
                                  mat_var(3)%matrix_diagonal_d,                    &
                                  mat_var(3)%packed_columns_d,                     &
                                  mat_var(3)%non_zero_columns,                     &
                                  mat_var(3)%row_index,                            &
                                  mat_var(3)%number,n3d,'hamiltonian',             &
                                  print_packed)
                 DEALLOCATE(triangle_hamiltonian_d)
             END IF
          END IF
      ELSE
!
!         Use Full Storage
!
          ALLOCATE(hamiltonian_d(n3d,n3d))
          IF(non_orth) THEN      
             ALLOCATE(overlap_d(n3d,n3d))
          END IF
          Call Compute_Eigenstates(hamiltonian_d, overlap_d) 
          IF(preconditioner == 'cholesky_decomposition') THEN
             ALLOCATE(upper_d(n3d*(n3d+1)/2), lower_d(n3d*(n3d+1)/2)) 
             write(iout,3)
             Call Full_to_Upper(overlap_d,upper_d)
             CALL Cholesky(upper_d,lower_d,print_parameter)
             IF (packed) THEN
!
!                Use packed storage from now on.
!
                 DEALLOCATE(upper_d, lower_d)
                 ALLOCATE(mat_var(3)%non_zero_columns(n3d),                        &
                          mat_var(3)%row_index(lenbuf),                            &
                          mat_var(3)%packed_columns_d(lenbuf),                     &
                          mat_var(3)%matrix_diagonal_d(n3d) )
                 Call pack_matrix(hamiltonian_d,                                   &
                                  mat_var(3)%matrix_diagonal_d,                    &
                                  mat_var(3)%packed_columns_d,                     &
                                  mat_var(3)%non_zero_columns,                     &
                                  mat_var(3)%row_index,                            &
                                  mat_var(3)%number,n3d,'hamiltonian',             &
                                  print_packed)
                 DEALLOCATE(hamiltonian_d)
                 ALLOCATE(mat_var(4)%non_zero_columns(n3d),                        &
                          mat_var(4)%row_index(lenbuf),                            &
                          mat_var(4)%packed_columns_d(lenbuf),                     &
                          mat_var(4)%matrix_diagonal_d(n3d) )
                 Call pack_matrix(overlap_d,                                       &
                                  mat_var(4)%matrix_diagonal_d,                    &
                                  mat_var(4)%packed_columns_d,                     &
                                  mat_var(4)%non_zero_columns,                     &
                                  mat_var(4)%row_index,                            &
                                  mat_var(4)%number,n3d,'overlap',                 &
                                  print_packed)
                 DEALLOCATE(overlap_d)
             END IF
          ELSE
             IF (packed) THEN
!
!                Use packed storage from now on.
!
                 ALLOCATE(mat_var(3)%non_zero_columns(n3d),                        &
                          mat_var(3)%row_index(lenbuf),                            &
                          mat_var(3)%packed_columns_d(lenbuf),                     &
                          mat_var(3)%matrix_diagonal_d(n3d) )
                 Call pack_matrix(hamiltonian_d,                                   &
                                  mat_var(3)%matrix_diagonal_d,                    &
                                  mat_var(3)%packed_columns_d,                     &
                                  mat_var(3)%non_zero_columns,                     &
                                  mat_var(3)%row_index,                            &
                                  mat_var(3)%number,n3d,'hamiltonian',             &
                                  print_packed)
                 DEALLOCATE(hamiltonian_d)
             END IF
          END IF
      END IF
  ELSE
      IF (triangle) THEN
          ALLOCATE(triangle_hamiltonian_z(n3d*(n3d+1)/2))
          IF(non_orth) THEN      
             ALLOCATE(triangle_overlap_z(n3d*(n3d+1)/2))
          END IF
          Call Compute_Eigenstates(triangle_hamiltonian_z, triangle_overlap_z) 
          IF(preconditioner == 'cholesky_decomposition') THEN
             ALLOCATE(upper_z(n3d*(n3d+1)/2), lower_z(n3d*(n3d+1)/2)) 
             write(iout,3)
             upper_z = triangle_overlap_z
             CALL Cholesky(upper_z,lower_z,print_parameter)
             IF (packed) THEN
!
!                Use packed storage from now on.
!
                 DEALLOCATE(upper_z, lower_z)
                 ALLOCATE(mat_var(3)%non_zero_columns(n3d),                        &
                          mat_var(3)%row_index(lenbuf),                            &
                          mat_var(3)%packed_columns_z(lenbuf),                     &
                          mat_var(3)%matrix_diagonal_z(n3d) )
                 Call h_pack_upper(triangle_hamiltonian_z,                         &
                                  mat_var(3)%matrix_diagonal_z,                    &
                                  mat_var(3)%packed_columns_z,                     &
                                  mat_var(3)%non_zero_columns,                     &
                                  mat_var(3)%row_index,                            &
                                  mat_var(3)%number,n3d,'hamiltonian',             &
                                  print_packed)
                 DEALLOCATE(triangle_hamiltonian_z)
                 ALLOCATE(mat_var(4)%non_zero_columns(n3d),                        &
                          mat_var(4)%row_index(lenbuf),                            &
                          mat_var(4)%packed_columns_z(lenbuf),                     &
                          mat_var(4)%matrix_diagonal_z(n3d) )
                 Call h_pack_upper(triangle_overlap_z,                             &
                                  mat_var(4)%matrix_diagonal_z,                    &
                                  mat_var(4)%packed_columns_z,                     &
                                  mat_var(4)%non_zero_columns,                     &
                                  mat_var(4)%row_index,                            &
                                  mat_var(4)%number,n3d,'overlap',                 &
                                  print_packed)
                 DEALLOCATE(triangle_overlap_z)
             END IF
          ELSE
             IF (packed) THEN
!
!                Use packed storage from now on.
!
                 ALLOCATE(mat_var(3)%non_zero_columns(n3d),                        &
                          mat_var(3)%row_index(lenbuf),                            &
                          mat_var(3)%packed_columns_z(lenbuf),                     &
                          mat_var(3)%matrix_diagonal_z(n3d) )
                 Call h_pack_upper(triangle_hamiltonian_z,                         &
                                  mat_var(3)%matrix_diagonal_z,                    &
                                  mat_var(3)%packed_columns_z,                     &
                                  mat_var(3)%non_zero_columns,                     &
                                  mat_var(3)%row_index,                            &
                                  mat_var(3)%number,n3d,'hamiltonian',             &
                                  print_packed)
                 DEALLOCATE(triangle_hamiltonian_z)
             END IF
          END IF
      ELSE
          ALLOCATE(hamiltonian_z(n3d,n3d))
          IF(non_orth) THEN      
             ALLOCATE(overlap_z(n3d,n3d))
          END IF
             Call Compute_Eigenstates(hamiltonian_z, overlap_z) 
             IF(preconditioner == 'cholesky_decomposition') THEN
                ALLOCATE(upper_z(n3d*(n3d+1)/2), lower_z(n3d*(n3d+1)/2)) 
                write(iout,3)
                Call Full_to_Upper(overlap_z,upper_z)
                CALL Cholesky(upper_z,lower_z,print_parameter)
                IF (packed) THEN
!
!                Use packed storage from now on.
!
                    DEALLOCATE(upper_z, lower_z)
                    ALLOCATE(mat_var(3)%non_zero_columns(n3d),                     &
                             mat_var(3)%row_index(lenbuf),                         &
                             mat_var(3)%packed_columns_z(lenbuf),                  &
                             mat_var(3)%matrix_diagonal_z(n3d) )
                    Call pack_matrix(hamiltonian_z,                                &
                                  mat_var(3)%matrix_diagonal_z,                    &
                                  mat_var(3)%packed_columns_z,                     &
                                  mat_var(3)%non_zero_columns,                     &
                                  mat_var(3)%row_index,                            &
                                  mat_var(3)%number,n3d,'hamiltonian',             &
                                  print_packed)
                    ALLOCATE(mat_var(4)%non_zero_columns(n3d),                     &
                             mat_var(4)%row_index(lenbuf),                         &
                             mat_var(4)%packed_columns_z(lenbuf),                  &
                             mat_var(4)%matrix_diagonal_z(n3d) )
                    Call pack_matrix(overlap_z,                                    &
                                  mat_var(4)%matrix_diagonal_z,                    &
                                  mat_var(4)%packed_columns_z,                     &
                                  mat_var(4)%non_zero_columns,                     &
                                  mat_var(4)%row_index,                            &
                                  mat_var(4)%number,n3d,'overlap',                 &
                                  print_packed)
                    DEALLOCATE(overlap_z)
                END IF
             ELSE
                IF (packed) THEN
!
!                Use packed storage from now on.
!
                    ALLOCATE(mat_var(3)%non_zero_columns(n3d),                     &
                             mat_var(3)%row_index(lenbuf),                         &
                             mat_var(3)%packed_columns_z(lenbuf),                  &
                             mat_var(3)%matrix_diagonal_z(n3d) )
                    Call pack_matrix(hamiltonian_z,                                &
                                  mat_var(3)%matrix_diagonal_z,                    &
                                  mat_var(3)%packed_columns_z,                     &
                                  mat_var(3)%non_zero_columns,                     &
                                  mat_var(3)%row_index,                            &
                                  mat_var(3)%number,n3d,'hamiltonian',             &
                                  print_packed)
                END IF
             END IF
      END IF
  END IF
  DEALLOCATE(rowlab,collab)
  IF(diagonalize_only) THEN
     return
  END IF
1 Format(/,20x,'Using Packed Forms for all Matrices')
2 FORMAT(/,10x,'Matrix Size = ',i8)
3 Format(/,20x,'Cholesky Decomposition of Overlap Matrix')
!***********************************************************************
  END SUBROUTINE Compute_Eigenstates_d
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
  END  MODULE Cholesky_Decomposition
!***********************************************************************
!***********************************************************************
