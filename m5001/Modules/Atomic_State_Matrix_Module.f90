!***********************************************************************
                           MODULE Atomic_State_Matrix_Module
!
                           USE accuracy
                           USE Data_Module,  ONLY : zero
                           USE Packed_Matrix_Module
                           USE Angular_Module
                           USE FEDVR_Shared, ONLY : file_loc
                           USE Packed_Matrix_Vector_Multiply_Module,                    &
                                             ONLY : u_l_tran_column_packed_matrix_u_r,  & 
                                                    u_l_tran_column_packed_matrix_u_r_on_vector 
!
                           IMPLICIT NONE
            CHARACTER(LEN=1600)                  :: data_card
            CHARACTER(LEN=80)                    :: pass
            REAL                                 :: elapsed_time            
!***********************************************************************
!***********************************************************************
                           Contains
!***********************************************************************
!***********************************************************************
!deck Pack_State_Matrices
!***begin prologue     Pack_State_Matrices
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Driver for generating packed matrices on a state-by-state,
!***                   basis as well as the diagonalization and transformation matrices 
!***                   of the atomic Hamiltonian.  These are required to perform the  
!***                   time-dependent propagation in a contracted atomic basis with
!***                   maximum efficiency.  At the conclusion of the subroutine one
!***                   has all of the Hamiltonian, overlap and dipole matrix elements
!***                   in packed form as well as the eigenvalues and eigenvectors for
!***                   each atomic L block.
!***
!***description        The input matrices for each angular momentum are are made up of 
!***                   channel-channel submatrices and channel-bound and bound-bound submatrices.  
!***                   Symbolically each looks like
!***                                   H_CC 
!***                                   H_BC H_BB
!***                   where only the triangle is stored when the matrices are symmetric,
!***                   such as the Hamiltonian and overlap.  The program removes unwanted splines 
!***                   associated with the channel matrices which, if left in, would lead to physically 
!***                   incorrect results.  Splines at the origin and at large distances are excised 
!***                   from the input and the matrices are stored in packed IOsys format.  For the dipole
!***                   matrices, which are in general rectangular, the same transformation is used to
!***                   remove the unwanted splines.  The matrices are then packed column-wise, with
!***                   all of the elements below a user-set threshold eliminated.
!***references
!***routines called
!***end prologue       Pack_State_Matrices
!
  SUBROUTINE Pack_State_Matrices
  IMPLICIT NONE
  REAL(idp)                             :: secnds
  REAL(idp)                             :: del_t
!
  DEALLOCATE(mat_var)
  IF(input_Matrices(24:43) == 'drake_channel_format') THEN
     time(1) = secnds(0.0) 
     write(iout,*)
     write(iout,*) '          *** Opening the State Matrix File ***'
     file_loc=File_Directory(4)(1:len_dir(4))//'/'//'state_matrix_file'
     Call pakstr(file_loc,len(1))
     Call IOsys('open state_matrices as new',0,0,0,file_loc(1:len(1)))
     ALLOCATE( labels(0:L_Max), state_matrix_size(0:L_Max),                            &
               state_mat(0:L_Max), dipole_mat(0:L_Max), matrix_size(0:L_Max) )
     Call Generate_Pointers
     write(iout,*)
     write(iout,*) '          *** Closing State Matrix File ***'
     Call IOsys('close state_matrices',0,0,0,' ')
     write(iout,*)
     write(iout,*) '          *** Opening the Packed Matrix File ***'
     file_loc=File_Directory(4)(1:len_dir(4))//'/'//'packed_matrix_file'
     Call pakstr(file_loc,len(1))
     Call IOsys('open packed_matrices as new',0,0,0,file_loc(1:len(1)))
     Call Reformat_State_Matrices
     write(iout,*)
     write(iout,*) '          *** Closing Packed Matrix File ***'
     Call IOsys('close packed_matrices',0,0,0,' ')
     DO L = 0, L_Max
        DEALLOCATE(labels(L)%state_quantum_numbers, state_mat(L)%matrix_pointer)
     END DO
     DEALLOCATE( labels, state_matrix_size, state_mat, matrix_size, dipole_mat )
     IF( reformat_input_matrix_only) THEN
         write(iout,*) '          *** Quitting After Reformat ***'
         time(2) = secnds(0.0)
         del_t = time(2) - time(1)
         WRITE(iout,1)
         WRITE(iout,2) del_t
         WRITE(iout,1)
         stop
     ELSE
         Call Input_Reformatted_State_Matrices
         time(2) = secnds(0.0)
         del_t = time(2) - time(1)
         WRITE(iout,1)
         WRITE(iout,3) del_t
         WRITE(iout,1)
         stop
     END IF
  ELSE IF(input_Matrices(24:45) == 'packed_in_iosys_format') THEN
         time(1) = secnds(0.0)
         stop
         Call Input_Reformatted_State_Matrices
         time(2) = secnds(0.0)
         del_t = time(2) - time(1)
         WRITE(iout,1)
         WRITE(iout,4) del_t
         WRITE(iout,1)
         stop
  ELSE
     Write(iout,5) input_matrices
     Call lnkerr('Bad Input Matrix Keyword')
  END IF
1 FORMAT('***********************************************'                               &
         '*************************')
2 FORMAT(/,10X,'Time to Reformat, Read, Write, Pack and Cholesky Decompose Matrices = ',f15.8)
3 FORMAT(/,10X,'Time to Reformat, Read, Write, Pack, Cholesky Decompose and Form Dipole' &
               ' Matrices = ',f15.8)
4 FORMAT(/,10X,'Time to Open and Read Packed Matrices and Form Dipole Matrices = ',f15.8)
5 FORMAT(/,10X,'Bad Input Matrix Keyword Used = ',a80,/,10x,                             &
               'Allowed Values are: using_angular_symmetry_drake_channel_format',/,10x,  &
               '                    using_angular_symmetry_packed_in_iosys_format')
END SUBROUTINE Pack_State_Matrices
!***********************************************************************
!***********************************************************************
!deck Input_Reformatted_State_Matrices
!***begin prologue     Input_Reformatted_State_Matrices
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Generate pointers to reformat matrices
!***
!***description        the input matrices are made up of channel-channel submatrices and
!***                   channel-bound and bound-bound submatrices.  symbolically it looks like
!***                                   H_CC 
!***                                   H_BC H_BB
!***                   only the triangle is stored where the matrices are symmetric such as the Hamiltonian
!***                   and overlap.  the program removes unwanted splines associated with the channel matrices
!***                   which, if left in, would lead to physically incorrect results.  splines at the origin and
!***                   at large distances are excised from the input and the matrices are stored in IOsys format.
!***references
!***routines called
!***end prologue       Input_Reformatted_State_Matrices
!
  SUBROUTINE Input_Reformatted_State_Matrices
  IMPLICIT NONE
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: IOSTAT
  INTEGER                               :: m_sub
  INTEGER                               :: count
  LOGICAL                               :: dollar
  REAL(idp), DIMENSION(:), ALLOCATABLE  :: vector_in
  REAL(idp), DIMENSION(:), ALLOCATABLE  :: vector_out
!
!
  write(iout,*)
  write(iout,*) '          *** Reopening File to Read State Matrices ***'
  file_loc=File_Directory(4)(1:len_dir(4))//'/'//'state_matrix_file'
  Call pakstr(file_loc,len(1))
  Call IOsys('open state_matrices as old',0,0,0,file_loc(1:len(1)))
  write(iout,*) '              *** Reading Needed Variables From State Matrices ***'
  Call iosys('read integer L_Max from state_matrices',1,L_Max,0,' ')
  ALLOCATE( state_matrix_size(0:L_Max), final_state_matrix_size(0:L_Max),                  &
            state_mat(0:L_Max), dipole_mat(0:L_Max) )
  Call IOsys('read integer state_matrix_size from state_matrices',L_max+1,                 &
              state_matrix_size, 0 , ' ' )       
  write(iout,*) '              *** Closing State Matrices File ***'
  Call IOsys('close state_matrices',0,0,0,' ')          
  write(iout,*)
  write(iout,*) '          *** Reopening Packed Matrix File ***'
  file_loc=File_Directory(4)(1:len_dir(4))//'/'//'packed_matrix_file'
  Call pakstr(file_loc,len(1))
  Call IOsys('open packed_matrices as old',0,0,0,file_loc(1:len(1)))
  write(iout,*) '              *** Reading Needed Variables From Packed Matrices ***'
  count = 0
  DO L = 0, L_Max
     file_key = 'L_'//itoc(L)
     ALLOCATE ( state_mat(L)%state_eigen_values(state_matrix_size(L)) )
     Call IOsys('read real "hamiltonian_eigenvalues_'//file_key//'" from '                 &
                //'packed_matrices', state_matrix_size(L),                                 &
                   state_mat(L)%state_eigen_values,0,' ')
     write(iout,*) '              *** Reading in Eigenvalues ***'
     final_state_matrix_size(L) = 0
     DO i = 1, state_matrix_size(L)
        IF ( state_mat(L)%state_eigen_values(i) < zero                                     &
                           .or.                                                            &
             abs(state_mat(L)%state_eigen_values(i)) < energy_cutoff ) THEN
             count = count + 1
             final_state_matrix_size(L) = final_state_matrix_size(L) + 1
        END IF
     END DO
     write(iout,*) ' L = ', L,' Number of Eigenvalues Retained = ',                        &
                                final_state_matrix_size(L)

  END DO
  Call IOsys('write integer final_state_matrix_size to packed_matrices',L_max+1,           &
              final_state_matrix_size, 0 , ' ' )       
  Call IOsys('write integer final_total_matrix_size to packed_matrices',1,                 &
              count, 0 , ' ' )       
  write(iout,*) ' Total Size of Final Matrix = ',count
  ALLOCATE( atomic_eigen_values(1:count)  )
  first = 0               
  DO L=0, L_Max
     first = first + 1
     last = first + final_state_matrix_size(L) - 1
     atomic_eigen_values(first:last) = state_mat(L)%state_eigen_values(1:final_state_matrix_size(L))
     first = last
     DEALLOCATE ( state_mat(L)%state_eigen_values )                
  END DO
  Call IOsys('write real "atomic_eigenvalues" to ' //'packed_matrices', count,             &
              atomic_eigen_values,0,' ')
  IF (Dipole_Matrices) THEN
      DO L=0,L_Max-1
!
!        Transform the dipole matrix from the original basis to the truncated set
!        of atomic eigenstates.
!
         local_time = secnds(0.0)
         ALLOCATE(state_mat(L)%state_eigen_vectors(state_matrix_size(L),                   &
                                                   final_state_Matrix_size(L)),            &
                  state_mat(L+1)%state_eigen_vectors(state_matrix_size(L+1),               &
                                                   final_state_Matrix_size(L+1)))
         write(iout,*) '              *** Reading in Eigenvectors ***'
         file_key = 'L_'//itoc(L+1)
         Call IOsys('read real "hamiltonian_eigenvectors_'//file_key//'" from '            &
                    //'packed_matrices',                                                   &
                       state_matrix_size(L+1)*final_state_matrix_size(L+1),                &
                       state_mat(L+1)%state_eigen_vectors,0,' ') 
         write(iout,*) ' Eigenvectors Read in for '//file_key
         file_key = 'L_'//itoc(L)
         Call IOsys('read real "hamiltonian_eigenvectors_'//file_key//'" from '            &
                    //'packed_matrices',                                                   &
                       state_matrix_size(L)*final_state_matrix_size(L),                    &
                       state_mat(L)%state_eigen_vectors,0,' ') 
         write(iout,*) ' Eigenvectors Read in for '//file_key
         Call IOsys('read integer "number_non_zero_dipole_matrix_'//file_key//             &
                    '_elements" from packed_matrices',1,                                   &
                      dipole_mat(L)%dipole_number,0,' ')
         ALLOCATE (dipole_mat(L)%dipole_non_zero_columns(state_matrix_size(L+1)),          &
                   dipole_mat(L)%dipole_row_index(dipole_mat(L)%dipole_number),            &
                   dipole_mat(L)%dipole_packed_columns(dipole_mat(L)%dipole_number) )                    
         Call Read_and_Write_Column_Packed_Matrices (                                      &
                                  dipole_mat(L)%dipole_packed_columns,                     &
                                  dipole_mat(L)%dipole_non_zero_columns,                   &
                                  dipole_mat(L)%dipole_row_index,                          &
                                  'read',                                                  &
                                  'dipole_matrix_'//file_key,                              &
                                  dipole_mat(L)%dipole_number)                       
         ALLOCATE( dipole_mat(L)%dipole_coupling_matrix( final_state_matrix_size(L),       &
                                                         final_state_matrix_size(L+1)))
         Call u_l_tran_column_packed_matrix_u_r (                                          &
                                                state_mat(L)%state_eigen_vectors,          &
                                                state_mat(L+1)%state_eigen_vectors,        &
                                                dipole_mat(L)%dipole_packed_columns,       &
                                                dipole_mat(L)%dipole_non_zero_columns,     &
                                                dipole_mat(L)%dipole_row_index,            &
                                                dipole_mat(L)%dipole_coupling_matrix )
         DEALLOCATE(state_mat(L)%state_eigen_vectors,                                      &
                    state_mat(L+1)%state_eigen_vectors)
!         ALLOCATE( vector_in(final_state_matrix_size(L+1)),                                &
!                   vector_out(final_state_matrix_size(L+1)) )
!         Call u_l_tran_column_packed_matrix_u_r_on_vector (                                &
!                                                state_mat(L)%state_eigen_vectors,          &
!                                                state_mat(L+1)%state_eigen_vectors,        &
!                                                vector_in,                                 &
!                                                vector_out,                                &
!                                                dipole_mat(L)%dipole_packed_columns,       &
!                                                dipole_mat(L)%dipole_non_zero_columns,     &
!                                                dipole_mat(L)%dipole_row_index )
         DEALLOCATE (dipole_mat(L)%dipole_non_zero_columns,                                &
                     dipole_mat(L)%dipole_row_index,                                       &
                     dipole_mat(L)%dipole_packed_columns )                    
         lenbuf = final_state_matrix_size(L) * final_state_matrix_size(L+1)
         ALLOCATE (dipole_mat(L)%dipole_non_zero_columns(final_state_matrix_size(L+1)),    &
                   dipole_mat(L)%dipole_row_index(lenbuf),                                 &
                   dipole_mat(L)%dipole_packed_columns(lenbuf) )                    
         drop_tol = drop_dipole
         Call Pack_Rectangular  (                                                           &
                                 dipole_mat(L)%dipole_coupling_matrix,                      &
                                 dipole_mat(L)%dipole_packed_columns,                       &
                                 dipole_mat(L)%dipole_non_zero_columns,                     &
                                 dipole_mat(L)%dipole_row_index,                            &
                                 dipole_mat(L)%dipole_number,                               &
                                 final_state_matrix_size(L),                                &
                                 final_state_matrix_size(L+1),                              &
                                'dipole' )
         write(iout,*)
         write(iout,*) '     Repacking trimmed dipole for L_1 = ', L, '  L_2 = ', L+1
         write(iout,*)
         write(iout,*) '     Left Dimension           = ' , final_state_matrix_size(L)      
         Write(iout,*) '     Right Dimension          = ' , final_state_matrix_size(L+1)    
         Write(iout,*) '     Non Zero Column Elements = ' , dipole_mat(L)%dipole_number     
         Write(iout,*) '     Writing                  = ' , 'dipole_matrix_'//file_key            
         Call Read_and_Write_Column_Packed_Matrices (                                      &
                                                    dipole_mat(L)%dipole_packed_columns,   &
                                                    dipole_mat(L)%dipole_non_zero_columns, &
                                                    dipole_mat(L)%dipole_row_index,        &
                                                    'write',                               &
                                                    'c_dipole_matrix_'//file_key,          &
                                                    dipole_mat(L)%dipole_number )
         DEALLOCATE( dipole_mat(L)%dipole_coupling_matrix,                                 &
                     dipole_mat(L)%dipole_non_zero_columns,                                &
                     dipole_mat(L)%dipole_row_index,                                       &
                     dipole_mat(L)%dipole_packed_columns )                    
         elapsed_time = secnds(0.0) - local_time
         write(iout,1) L, L+1, elapsed_time
     END DO
  END IF
  write(iout,*)
  write(iout,*) '          *** Closing Packed Matrices ***'
  Call IOsys('close packed_matrices',0,0,0,' ')          
!
1 Format(/,10x,'Time to Transform and Re-pack the Dipole Matrices in the Atomic Basis'     &
         /,10x,'for L = ',i3,2x,', L+1 =',i3,                                              &
         /,10x,'      = ',f15.8)
END SUBROUTINE Input_Reformatted_State_Matrices
!***********************************************************************
!***********************************************************************
!deck Generate_Pointers
!***begin prologue     Generate_Pointers
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Generate pointers to reformat matrices
!***
!***description        the input matrices are made up of channel-channel submatrices and
!***                   channel-bound and bound-bound submatrices.  symbolically it looks like
!***                                   H_CC 
!***                                   H_BC H_BB
!***                   only the triangle is stored where the matrices are symmetric such as the Hamiltonian
!***                   and overlap.  the program removes unwanted splines associated with the channel matrices
!***                   which, if left in, would lead to physically incorrect results.  splines at the origin and
!***                   at large distances are excised from the input and the matrices are stored in IOsys format.
!***references
!***routines called
!***end prologue       Generate_Pointers
!
  SUBROUTINE Generate_Pointers
  IMPLICIT NONE
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: IOSTAT
  INTEGER                               :: m_sub
  LOGICAL                               :: dollar
  INTEGER                               :: count
!

!
  Call iosys('write integer L_Max to state_matrices',1,L_Max,0,' ')
!
! First just create labels and a pointer to all of the matrix indices that will be       
! retained in the problem.
!                         Do not read the matrices at this point.
!
  local_time = secnds(0.0)
  state_matrix_size(0:L_Max) = 0
  n3d = 0
  DO L = 0, L_Max
     L_1=itoc(L)
     IF ( dollar('$state_data_l_'//L_1,data_card,pass,inp) ) THEN
          matrix_file=chrkey(data_card,'h_matrix_file_name_l='//L_1,                               &
                                       'bsr_mat.'//L_1,' ')
          Call pakstr(matrix_file,len(1))
          parity=chrkey(data_card,'parity','even',' ')
          file_loc=File_Directory(4)(1:len_dir(4))//'/'//matrix_file(1:len(1))
          Call pakstr(file_loc,len(2))
          OPEN(UNIT=50,FILE=file_loc(1:len(2)),ACCESS='sequential',FORM='unformatted', IOSTAT=IOSTAT,STATUS='old')
          READ(50) number_of_splines, number_of_channels, number_of_correlation_terms
          CLOSE(50)
          matrix_size(L) = number_of_splines * number_of_channels + number_of_correlation_terms
          n3d = n3d + matrix_size(L)
          ALLOCATE(labels(L)%state_quantum_numbers(0:number_of_channels+2) )
          labels(L)%state_quantum_numbers(0) = 0
          IF (parity == 'odd') THEN
              labels(L)%state_quantum_numbers(0) = 1
          END IF
          labels(L)%state_quantum_numbers(1) = number_of_channels
          labels(L)%state_quantum_numbers(2) = number_of_correlation_terms
          Call intarr(data_card,'channel_l_values',labels(L)%state_quantum_numbers(3),             &
                      number_of_channels,' ')
          Call IOsys('write integer "channel_labels_state = '//L_1//'" to state_matrices',         &
                      number_of_channels+3,labels(L)%state_quantum_numbers,0, ' ' )             
          ALLOCATE( state_mat(L)%matrix_pointer(matrix_size(L) ))
          first_matrix_index = 0
          last_matrix_index = 0
          count = 0
          write(iout,1) L
          DO i = 1, number_of_channels       
             first_matrix_index = last_matrix_index + 1
             last_matrix_index  = last_matrix_index + number_of_splines
!
!            The pointer to the first index of the original matrix to be used 
!            in the channel block i
!
             remove = labels(L)%state_quantum_numbers(2 + i) + 1
             remove = min(spline_order-1,remove)
             first_new_matrix_index = first_matrix_index + remove
!         
!            The pointer to the last index of the original matrix to be used 
!            in the channel block i
!
             IF(remove_last_spline) THEN
                last_new_matrix_index = last_matrix_index - 1
             END IF
             IF(remove_next_to_last_spline) THEN
                last_new_matrix_index = last_new_matrix_index  - 1
             END IF
!            Calculate the pointers
!
             DO j = first_new_matrix_index, last_new_matrix_index
                count = count + 1        
                state_mat(L)%matrix_pointer(count) = j
             END DO
             m_sub = last_new_matrix_index - first_new_matrix_index + 1
             write(iout,2) labels(L)%state_quantum_numbers( 2+i),                                  &
                           m_sub, first_new_matrix_index, last_new_matrix_index
!
!            Accumulate size of the state.
!
             state_matrix_size(L) = state_matrix_size(L) + last_new_matrix_index                   &
                                                           -                                       &
                                    first_new_matrix_index  + 1
!
          END DO
!
!
!        The correlation terms if present
         IF ( number_of_correlation_terms /= 0 ) THEN
              Write(iout,5)
              DO i = last_matrix_index + 1 , last_matrix_index + number_of_correlation_terms
                 count = count + 1
                 state_mat(L)%matrix_pointer(count) = i
              END DO
              Write(iout,6) number_of_correlation_terms, last_matrix_index + 1,                    &
                            last_matrix_index + number_of_correlation_terms
         END IF
!
!        The full size of this L block after removal of the unwanted splines and adding in the 
!        correlation terms.
!
         state_matrix_size(L) = state_matrix_size(L) + number_of_correlation_terms
         Write(iout,3) matrix_file(1:len(1)), number_of_channels,                                  &
                       number_of_correlation_terms, state_matrix_size(L),                          &
                       parity
         IF(write_pointers) THEN
            write(iout,4) state_mat(L)%matrix_pointer(1:count)
         END IF 
         Call IOsys('write integer "pointers_state = '//L_1//'" to state_matrices',                &
                     count,state_mat(L)%matrix_pointer, 0, ' ' )                 
     END IF
  END DO
  elapsed_time = secnds(0.0) - local_time
  write(iout,7) elapsed_time
  Call IOsys('write integer matrix_size to state_matrices', L_Max + 1, matrix_size, 0 , ' ' )  
  Call IOsys('write integer state_matrix_size to state_matrices', L_Max + 1,                       &
              state_matrix_size, 0 , ' ' )                 
1 Format(/,30x,'Matrix Information for L = ',i3,/,10x,'Channel Angular Momentum',                  &
            5x,'Submatrix Size',5x,'First Index',5x,'Last Index')
2 Format(20x,i3,20x,i5,12x,i5,10x,i5)
3 Format(/10x,'Opening Hamiltonian File    = ',a24,/,15x,                                          &
              'Number of Channels          = ',i4,/,15x,                                           &
              'Number of Correlation Terms = ',i4,/,15x,                                           &
              'State Matrix Size           = ',i4,/,15x,                                           &
              'Parity                      = ',a4 )
4 Format(/,15x,5X,'Pointer = ',(/,20i5))
5 Format(/,10x,'Number of Correlation Terms', 21x,'First Index',5x,'Last Index')
6 Format(20x,i5,35x,i5,10x,i5)
7 Format(/,10x,'Time to Generate the Pointers = ',f15.8)
END SUBROUTINE Generate_Pointers
!***********************************************************************
!***********************************************************************
!deck Reformat_State_Matrices
!***begin prologue     Reformat_State_Matrices
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Read input from B-spline program and reformat matrices
!***
!***description        The input matrices are made up of channel-channel submatrices and
!***                   channel-bound and bound-bound submatrices. Symbolically 
!***                   it looks like
!***                                   H_CC 
!***                                   H_BC H_BB
!***                   where only the triangle is stored when the matrices are symmetric 
!***                   such as the Hamiltonian and overlap.  The input matrices do not
!**                    have the zeros or near zeros removed.   The program removes 
!***                   unwanted splines associated with the channel matrices
!***                   which, if left in, would lead to physically incorrect results.  
!***                   Splines at the origin and at large distances are excised 
!***                   from the input and the output matrices are stored in packed IOsys 
!***                   format using a column packing format where zeros are removed.  
!***                   There is one matrix for each L value in the case of the overlap and Hamiltonian.  For the
!***                   transition dipole matrix, we store the matrices coupling a given L
!***                   to its (L+1) and (L-1) partners.    
!***references
!***routines called
!***end prologue          Reformat_State_Matrices
!
  SUBROUTINE Reformat_State_Matrices
  IMPLICIT NONE
  TYPE(REAL_MATRICES)                   :: mat_d
  CHARACTER (LEN=8)                     :: itoc
  INTEGER                               :: lenth
  INTEGER                               :: IOSTAT
  INTEGER                               :: lwork
!
!
! Now go trim the actual overlap and Hamiltonian matrices
!
  DO L=0, L_Max
     L_1 = itoc(L)
     matrix_file='bsr_mat.'//L_1
     L_Prime = L + 1
     L_2 = itoc(L_Prime)
     Call pakstr(matrix_file,len(1))
     file_loc=File_Directory(4)(1:len_dir(4))//'/'//matrix_file(1:len(1))
     Call pakstr(file_loc,len(2))
     write(iout,*)
     write(iout,*) '     Open the Hamiltonian Input File = ', matrix_file(1:len(1))
     OPEN(UNIT=50,FILE=file_loc(1:len(2)),ACCESS='sequential',FORM='unformatted',IOSTAT=IOSTAT,STATUS='old')
     READ(50) 
     write(iout,*)
     write(iout,*) '     Trim and Pack Overlap and Hamiltonian Matrices'
     file_key = 'L_'//L_1
     tri_size = matrix_size(L) * ( matrix_size(L) + 1) / 2
     state_tri_size = state_matrix_size(L) * ( state_matrix_size(L) + 1) / 2
     lenbuf = state_tri_size - state_matrix_size(L)
     ALLOCATE (state_mat(L)%state_s_matrix(state_tri_size) )
     ALLOCATE (state_mat(L)%s_non_zero_columns(state_matrix_size(L)),                              &
               state_mat(L)%s_row_index(lenbuf),                                                   &
               state_mat(L)%s_packed_columns(lenbuf),                                              &
               state_mat(L)%s_diagonal(state_matrix_size(L) ) )
     local_time = secnds(0.0)
     drop_tol = drop_overlap
     Call Trim_and_Pack_Triangles(state_mat(L)%matrix_pointer,                                     &
                                  state_mat(L)%state_s_matrix,                                     &
                                  state_mat(L)%s_non_zero_columns,                                 &
                                  state_mat(L)%s_row_index,                                        &
                                  state_mat(L)%s_packed_columns,                                   &
                                  state_mat(L)%s_diagonal,                                         &
                                  state_mat(L)%s_number,                                           &
                                  matrix_size(L),                                                  &
                                  state_matrix_size(L),                                            &
                                  'overlap')
     IF (diagonalize_matrix == .false.) THEN
         DEALLOCATE ( state_mat(L)%state_s_matrix )
     END IF
     DEALLOCATE (state_mat(L)%s_non_zero_columns,                                                  &
                 state_mat(L)%s_row_index,                                                         &
                 state_mat(L)%s_packed_columns,                                                    &
                 state_mat(L)%s_diagonal )
!
     ALLOCATE ( state_mat(L)%state_h_matrix(state_tri_size) )
     ALLOCATE ( state_mat(L)%h_non_zero_columns(state_matrix_size(L)),                             &
                state_mat(L)%h_row_index(lenbuf),                                                  &
                state_mat(L)%h_packed_columns(lenbuf),                                             &
                state_mat(L)%h_diagonal(state_matrix_size(L) ) )
     drop_tol = drop_hamiltonian
     Call Trim_and_Pack_Triangles(state_mat(L)%matrix_pointer,                                     &
                                  state_mat(L)%state_h_matrix,                                     &
                                  state_mat(L)%h_non_zero_columns,                                 &
                                  state_mat(L)%h_row_index,                                        &
                                  state_mat(L)%h_packed_columns,                                   &
                                  state_mat(L)%h_diagonal,                                         &
                                  state_mat(L)%h_number,                                           &
                                  matrix_size(L),                                                  &
                                  state_matrix_size(L),                                            &
                                  'hamiltonian')
     elapsed_time = secnds(0.0) - local_time
     write(iout,3) L, state_matrix_size(L), elapsed_time
     IF (diagonalize_matrix == .false.) THEN
         DEALLOCATE ( state_mat(L)%state_h_matrix )
     END IF
     DEALLOCATE ( state_mat(L)%h_non_zero_columns,                                                 &
                  state_mat(L)%h_row_index,                                                        &
                  state_mat(L)%h_packed_columns,                                                   &
                  state_mat(L)%h_diagonal )
     REWIND(50)
     write(iout,*)
     write(iout,*) '     Close the Hamiltonian File'
     CLOSE(UNIT=50)
     IF (diagonalize_matrix == .true.) THEN
         local_time = secnds(0.0)
         Call Diagonalize_Packed_Matrices(  mat_d,                                                 &
                                            'input',                                               &
                                             state_matrix_size(L),                                 &
                                             .true.,                                               &
                                             .true.,                                               &
                                             file_key,                                             &
                                             state_mat(L)%state_h_matrix,                          &
                                             state_mat(L)%state_s_matrix )
         elapsed_time = secnds(0.0) - local_time
         write(iout,4) L, state_matrix_size(L), elapsed_time
         DEALLOCATE ( state_mat(L)%state_h_matrix, state_mat(L)%state_s_matrix )
     END IF
  END DO
  IF (Dipole_Matrices) THEN
      DO L=0, L_Max - 1
         L_1 = itoc(L)
         L_Prime = L + 1
         L_2 = itoc(L_Prime)
         matrix_file='dd_mat.'//L_1
         Call pakstr(matrix_file,len(1))
         matrix_file=matrix_file(1:len(1))//'.'//L_2
         Call pakstr(matrix_file,len(1))
         file_loc=File_Directory(4)(1:len_dir(4))//'/'//matrix_file(1:len(1))
         Call pakstr(file_loc,len(2))
         write(iout,*)
         write(iout,*) '     Open the Dipole Input File = ', matrix_file(1:len(1))
         OPEN(UNIT=50,FILE=file_loc(1:len(2)),ACCESS='sequential',FORM='unformatted',IOSTAT=IOSTAT,STATUS='old')
         local_time = secnds(0.0)
         READ(50) 
         write(iout,*)
         write(iout,*) '     Trim and Pack Dipole File'
         file_key = 'L_'//itoc(L)
         lenbuf = state_matrix_size(L) * state_matrix_size(L+1)
         ALLOCATE( dipole_mat(L)%dipole_coupling_matrix( state_matrix_size(L),                     &
                                                         state_matrix_size(L+1)))
         ALLOCATE (dipole_mat(L)%dipole_non_zero_columns(state_matrix_size(L+1)),                  &
                   dipole_mat(L)%dipole_row_index(lenbuf),                                         &
                   dipole_mat(L)%dipole_packed_columns(lenbuf) )
         drop_tol = drop_dipole
         Call Trim_and_Pack_Rectangles(state_mat(L)%matrix_pointer,                                &
                                       state_mat(L+1)%matrix_pointer,                              &
                                       dipole_mat(L)%dipole_coupling_matrix,                       &
                                       dipole_mat(L)%dipole_non_zero_columns,                      &
                                       dipole_mat(L)%dipole_row_index,                             &
                                       dipole_mat(L)%dipole_packed_columns,                        &
                                       dipole_mat(L)%dipole_number,                                &
                                       matrix_size(L),matrix_size(L+1),                            &
                                       state_matrix_size(L),state_matrix_size(L+1),                &
                                       'dipole')
         DEALLOCATE(dipole_mat(L)%dipole_coupling_matrix )
         DEALLOCATE (dipole_mat(L)%dipole_non_zero_columns,                                        &
                     dipole_mat(L)%dipole_row_index,                                               &
                     dipole_mat(L)%dipole_packed_columns )      
         elapsed_time = secnds(0.0) - local_time
         write(iout,5) L, L + 1, state_matrix_size(L), state_matrix_size(L+1), elapsed_time
!**************************************************************************************************
!
!                             Test Code
!
!         ALLOCATE( dipole_mat(L)%dipole_coupling_matrix( state_matrix_size(L),                     &
!                                                         state_matrix_size(L+1)))
!         ALLOCATE (dipole_mat(L)%dipole_non_zero_columns(state_matrix_size(L+1)),                  &
!                   dipole_mat(L)%dipole_row_index(dipole_mat(L)%dipole_number),                    &
!                   dipole_mat(L)%dipole_packed_columns(dipole_mat(L)%dipole_number) )                    
!         Call Read_and_Write_Column_Packed_Matrices (                                              &
!                                                 dipole_mat(L)%dipole_packed_columns,              &
!                                                 dipole_mat(L)%dipole_non_zero_columns,            &
!                                                 dipole_mat(L)%dipole_row_index,                   &
!                                                 'read',                                           &
!                                                 'dipole_matrix_'//file_key,                       &
!                                                  dipole_mat(L)%dipole_number )                         
!         Call Fill_General_Matrix_from_Column_Buffers (                                            &
!                                                       dipole_mat(L)%dipole_coupling_matrix,       &
!                                                       dipole_mat(L)%dipole_packed_columns,        &
!                                                       dipole_mat(L)%dipole_non_zero_columns,      &
!                                                       dipole_mat(L)%dipole_row_index)
!         local_title= 'Printing Matrix from Column_Buffers'
!         call prntfmn(local_title,dipole_mat(L)%dipole_coupling_matrix,                            &
!                            state_matrix_size(L),                                                  &
!                            state_matrix_size(L+1),                                                &
!                            state_matrix_size(L),                                                  &
!                            state_matrix_size(L+1),                                                &
!                            iout,'e')
!         DEALLOCATE(dipole_mat(L)%dipole_coupling_matrix )
!         DEALLOCATE (dipole_mat(L)%dipole_non_zero_columns,                                        &
!                     dipole_mat(L)%dipole_row_index,                                               &
!                     dipole_mat(L)%dipole_packed_columns )                               
!**************************************************************************************************
         REWIND(50)
         write(iout,*)
         write(iout,*) '     Close the Dipole File'
         CLOSE(UNIT=50)
      END DO
  END IF
!
1  Format(/,25x,'Reformatting B Spline Matrices: L_Max    = ',i3,/,25x,                            &
               '                       Number of Splines  = ',i4,/,25x,                            &
               '                       Spline Order       = ',i2,/25x,                             &
               '                       Remove_Last_Spline = ',l1,/25x,                             &
               '                       Remove_Next_To_Last_Spline = ',l1)
2  Format(/25x,'L = ',i2,/,15x,'Opening Hamiltonian File  = ',a24,/,15x,                           &
                               'Number of Channels          = ',i4,/,15x,                          &
                               'Number of Correlation Terms = ',i4)
3 Format(/,10x,'Time to Trim Sub-Matrix for L = ',i3,2x,'Size = ',i6,                              &
         /,10x,'                              = ',f15.8)
4 Format(/,10x,'Time to Diagonalize Sub-Matrix L = ', i3,2x'Size = ', i6,                          &
         /,10x,'                                 = ',f15.8)
5 Format(/,10x,'Time to Trim Dipole Sub-Matrix L = ', i3,2x,', L+1 = ',i3,                         &
         /,10x,'Left Dimension = ',i6,2x,'Right Dimension = ',i6,                                  &
         /,10x,'               = ',f15.8)
END SUBROUTINE Reformat_State_Matrices
!***********************************************************************
!***********************************************************************
!deck Trim_and_Pack_Triangles
!***begin prologue     Trim_and_Pack_Triangles
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            reformat input hamiltonian matrix
!***
!***references
!***routines called
!***end prologue       Trim_and_Pack_Triangles
!
  SUBROUTINE Trim_and_Pack_Triangles(pointer,matrix,non_zero_columns,        &
                                     row_index,packed_columns,diagonal,      &
                                     number,initial_dim,final_dim,type)
  IMPLICIT NONE
  INTEGER, DIMENSION(:)                        :: pointer
  REAL(idp), DIMENSION(:)                      :: matrix
  INTEGER, DIMENSION(:)                        :: non_zero_columns
  INTEGER, DIMENSION(:)                        :: row_index
  REAL(idp), DIMENSION(:)                      :: packed_columns
  REAL(idp), DIMENSION(:)                      :: diagonal
  INTEGER                                      :: number
  INTEGER                                      :: initial_dim
  INTEGER                                      :: final_dim
  CHARACTER(LEN=*)                             :: type
!
  ALLOCATE(scratch_tri(tri_size) )
  Call pakstr(type,len(1))
  write(iout,*)
  write(iout,*) '     Remove Unwanted Splines for '//type(1:len(1))//' L = ', L
  Call Reformat_Hamiltonian_Matrix(                                             &
                                   scratch_tri,                                 &
                                   matrix,                                      &
                                   pointer,                                     &
                                   initial_dim,                                 &
                                   final_dim )
  DEALLOCATE( scratch_tri )
  Call h_pack_upper (                                                           &
                     matrix,                                                    &
                     diagonal,                                                  &
                     packed_columns,                                            & 
                     non_zero_columns,                                          &
                     row_index,                                                 &
                     number,                                                    &
                     final_dim,                                                 &
                     type )
  write(iout,*)
  write(iout,*) '     Packing '//type(1:len(1))//' for L = ', L
  write(iout,*)
  write(iout,*) '     Initial Dimension          = ' , initial_dim
  Write(iout,*) '     Final_Dimension            = ' , final_dim
  Write(iout,*) '     Non Zero Column Elements   = ' , number
  Write(iout,*) '     Writing                    = '//type//'_matrix_'//file_key            
  Call Read_and_Write_Column_Packed_Matrices (                                  &
                                   packed_columns,                              &
                                   non_zero_columns,                            &
                                   row_index,                                   &
                                   'write',                                     &  
                                   type//'_matrix_'//file_key,                  &
                                   number,                                      &
                                   diagonal )

  END SUBROUTINE Trim_and_Pack_Triangles
!***********************************************************************
!***********************************************************************
!deck Trim_and_Pack_Rectangles
!***begin prologue     Trim_and_Pack_Rectangles
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            reformat input hamiltonian matrix
!***
!***references
!***routines called
!***end prologue       Trim_and_Pack_Rectangles
!
  SUBROUTINE Trim_and_Pack_Rectangles(pointer_l,pointer_r,matrix,               &
                                     non_zero_columns, row_index,               &
                                     packed_columns,number,                     &
                                     initial_dim_l,initial_dim_r,               &
                                     final_dim_l,final_dim_r,type)
  INTEGER, DIMENSION(:)                        :: pointer_l
  INTEGER, DIMENSION(:)                        :: pointer_r
  REAL(idp), DIMENSION(:,:)                    :: matrix
  INTEGER, DIMENSION(:)                        :: non_zero_columns
  INTEGER, DIMENSION(:)                        :: row_index
  REAL(idp), DIMENSION(:)                      :: packed_columns
  INTEGER                                      :: number
  INTEGER                                      :: initial_dim_l
  INTEGER                                      :: initial_dim_r
  INTEGER                                      :: final_dim_l
  INTEGER                                      :: final_dim_r
  CHARACTER(LEN=*)                             :: type

!
  ALLOCATE( scratch_matrix( initial_dim_l, initial_dim_r ) )
  Call pakstr(type,len(1))
  write(iout,*)
  write(iout,*) '     Remove Unwanted Splines for '//type(1:len(1))//' L_1 = ',   &
                      L, '  L_2 = ', L+1
  Call Reformat_Dipole_Matrix(                                                    &
                              scratch_matrix,                                     &
                              matrix,                                             &
                              pointer_l,                                          &
                              pointer_r,                                          &
                              initial_dim_l,                                      &
                              initial_dim_r,                                      &
                              final_dim_l,                                        &
                              final_dim_r )
!
!
  Call Pack_Rectangular(matrix,                                                   &
                        packed_columns,                                           &
                        non_zero_columns,                                         &
                        row_index,                                                &
                        number,                                                   &
                        final_dim_l,                                              &
                        final_dim_r,                                              &
                        type)
  write(iout,*)
  write(iout,*) '     Packing '//type(1:len(1))//' for L_1 = ', L, '  L_2 = ', L+1
  write(iout,*)
  write(iout,*) '     Initial Left Dimension   = ' , initial_dim_l
  write(iout,*) '     Final Left Dimension     = ' , final_dim_l
  Write(iout,*) '     Initial Right Dimension  = ' , initial_dim_r
  Write(iout,*) '     Final Right Dimension    = ' , final_dim_r
  Write(iout,*) '     Non Zero Column Elements = ' , number
  Write(iout,*) '     Writing                  = '//type//'_matrix_'//file_key            
  Call Read_and_Write_Column_Packed_Matrices (                                    &
                           packed_columns,                                        &
                           non_zero_columns,                                      &
                           row_index,                                             &
                           'write',                                               &
                           type//'_matrix_'//file_key,                            &
                           number )

END SUBROUTINE Trim_and_Pack_Rectangles
!***********************************************************************
!***********************************************************************
!deck Reformat_Hamiltonian_Matrix
!***begin prologue     Reformat_Hamiltonian_Matrix
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            reformat input hamiltonian matrix
!***
!***references
!***routines called
!***end prologue       Reformat_Hamiltonian_Matrix
!
  SUBROUTINE Reformat_Hamiltonian_Matrix(input_matrix,output_matrix,           &
                                         matrix_pointer,size_in,size_out)
  IMPLICIT NONE
  INTEGER                                      :: size_in
  INTEGER                                      :: size_out
  INTEGER, DIMENSION(:)                        :: matrix_pointer
  REAL(idp), DIMENSION(:)                      :: input_matrix
  REAL(idp), DIMENSION(:)                      :: output_matrix
  INTEGER                                      :: i
  INTEGER                                      :: j
  INTEGER                                      :: ii
  INTEGER                                      :: jj
  INTEGER                                      :: i_tri
  INTEGER                                      :: ind
  INTEGER                                      :: count
  INTEGER                                      :: start
  INTEGER                                      :: finish
!
  DO i = 1 , size_in
     start = i * ( i - 1) /2 + 1
     finish = start + i - 1 
     READ(50) input_matrix(start : finish )
  END DO
!
  IF (print_cc ) THEN
      local_title='Input Triangle of Matrix'
      write(iout,1) local_title
      DO i=1, size_in
         start = i * ( i - 1) /2 + 1
         finish = start + i - 1 
         write(iout,2) i
         write(iout,3) input_matrix(start : finish)
     END DO       
  END IF
!
! Remove Unwanted Matrix Elements and Fill the state matrices
!
  count = 0 
  DO i = 1, size_out
     ii = matrix_pointer(i)
     i_tri = ii * ( ii - 1) /2
     DO j = 1, i
        count = count + 1
        ind = i_tri + matrix_pointer(j)
        output_matrix(count) = input_matrix(ind)                   
     END DO
  END DO
  IF (print_cc ) THEN
      local_title='Reformatted Output Triangle of Matrix'
      write(iout,1) local_title
      DO i=1, size_out
         start = i * ( i - 1) /2 + 1
         finish = start + i - 1 
         write(iout,2) i
         write(iout,3) output_matrix(start : finish)
     END DO       
  END IF
1 Format(a80)
2 Format(/,5x,'Row = ',i4)
3 Format( (15x,5e15.8) )
END SUBROUTINE Reformat_Hamiltonian_Matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!***********************************************************************
!deck Reformat_Dipole_Matrix
!***begin prologue     Reformat_Dipole_Matrix
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            reformat input matrices
!***
!***references
!***routines called
!***end prologue       Reformat_Dipole_Matrix
!
  SUBROUTINE Reformat_Dipole_Matrix(input_matrix,output_matrix,          &
                                    matrix_pointer_left,                 &
                                    matrix_pointer_right,                &
                                    size_left,size_right,                &
                                    new_size_left,new_size_right)
  IMPLICIT NONE
  INTEGER                                      :: size_left
  INTEGER                                      :: size_right
  INTEGER                                      :: new_size_left
  INTEGER                                      :: new_size_right
  REAL(idp), DIMENSION(:,:)                    :: input_matrix
  REAL(idp), DIMENSION(:,:)                    :: output_matrix
  REAL(idp)                                    :: F_3J
  INTEGER, DIMENSION(:)                        :: matrix_pointer_left
  INTEGER, DIMENSION(:)                        :: matrix_pointer_right
  INTEGER                                      :: i
  INTEGER                                      :: j
  INTEGER                                      :: ii
  INTEGER                                      :: jj
  INTEGER                                      :: prefac
!
  DO i=1,size_right
     READ(50) (input_matrix(j,i), j = 1, size_left)
  END DO
!
  IF (print_cc ) THEN
      local_title='Input Matrix'
      call prntfmn(local_title,input_matrix,size_left,size_right,        &
                   size_left,size_right,iout,'e')
  END IF
!
! Remove Unwanted Matrix Elements and Fill the state matrices
!
  DO i = 1, new_size_right
     ii = matrix_pointer_right(i)
     DO j = 1, new_size_left 
        jj = matrix_pointer_left(j)        
        output_matrix(j,i) = input_matrix(jj,ii)                   
     END DO
  END DO
!
!             Go from reduced dipole matrix element to proper dipole matrix
!             elements using 3_J coefficients.  Here for linear polarization.  
!
  prefac = 1
  IF ( (L - 2 * L/2 ) == 1 ) THEN
        prefac = - 1
  END IF
  output_matrix(:,:) =  Three_J_Coef(L,0,1,0,L+1,0,.false.)             &
                            *                                           &
                        output_matrix(:,:)                              &
                            *                                           &
                          prefac
!
  IF (print_cc ) THEN
      local_title='Output Reformatted Matrix'
      call prntfmn(local_title,output_matrix,new_size_left,             &
                   new_size_right,new_size_left,new_size_right,iout,'e')
  END IF
!

END SUBROUTINE Reformat_Dipole_Matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!***********************************************************************
  END  MODULE Atomic_State_Matrix_Module
!***********************************************************************
!***********************************************************************
