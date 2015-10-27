!***********************************************************************
                           MODULE Packed_Matrix_Module
!
                           USE Atomic_Matrices
                           USE Global_Time_Propagation_Module
                           USE Pack_Hamiltonian_Module
                           USE Preconditioner_Module
!
                           IMPLICIT NONE
            CHARACTER(LEN=80)                    :: local_title
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                               Contains
!***********************************************************************
!***********************************************************************
!deck Packed_Matrix
!***begin prologue     Packed_Matrix
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Hamiltonian input and manipulation.
!***
!***description        This routine can reformat an input Hamiltonian
!                      into a packed form with zeros removed or read and
!                      store already packed matrices for use by the main
!                      routines.  Note that the packed matrices are stored in 
!                      IOsys format for easy reading and writing.
!***references
!***routines called
!***end prologue       Packed_Matrix
!
  SUBROUTINE Packed_Matrix(ham_type,file_directory)
  IMPLICIT NONE
  CHARACTER(LEN=*)                      :: ham_type
  CHARACTER(LEN=*)                      :: file_directory
  INTEGER                               :: lenth
  INTEGER                               :: len
  INTEGER                               :: len_1
  INTEGER                               :: IOSTAT
  REAL*4                                :: secnds
  REAL*4                                :: del_t
!
!
  IF (input_matrices(11:22) == 'iosys_format') THEN
      write(iout,*)
      write(iout,*) '          *** Reopening File to Read Packed Matrices ***'
      time(1) = secnds(0.0)
      len=lenth(file_directory)
      len_1=lenth(packed_file_name)
      Call IOsys('open packed_matrices as old',0,0,0,                            &
                  file_directory(1:len)//'/'//packed_file_name(1:len_1))
      Call Input_Packed_IOsys_Matrix(ham_type,file_directory)
      write(iout,*) '***Closing File Packed Matrices ***'
      Call IOsys('close packed_matrices',0,0,0,' ')
      time(2) = secnds(0.0)
      del_t = time(2) - time(1)
      WRITE(iout,1)
      WRITE(iout,3) del_t
      WRITE(iout,1)
  ELSE
      write(iout,*)
      write(iout,*) '          Reformatting Packed Matrices to IOsys Format'
!
!        In this branch of the if, we come in with a packed matrix from the Drake codes, 
!        and go out with a packed matrix in IOsys format.  The packed matrix from IOsys 
!        is done in the standard or triangular manner.  In standard manner an array of 
!        two indices and an array of elements that correspond to those 
!        indices is output along with the value of the number of non zero elements.  
!        Only the lower triangle is packed.  The other alternative is a packing of the 
!        triangle of the matrix.  This column oriented packing requires more extensive
!        logic but is needed to efficiently solve the linear equations when a 
!        non-orthogonal basis is used.
!
      time(1) = secnds(0.0) 
      len=lenth(file_directory)
      len_1=lenth(ham_file)
      write(iout,*) '          *** Opening and Processing Input Matrix File ***'
      OPEN(UNIT=50,FILE=file_directory(1:len)//'/'//ham_file(1:len_1),               &
           ACCESS='sequential', FORM='unformatted',IOSTAT=IOSTAT,STATUS='old')
      READ(50)
      write(iout,*) '          *** Opening File to Write Packed Matrices ***'
      len_1=lenth(packed_file_name)
      Call IOsys('open packed_matrices as new',0,0,0,                                &
                  file_directory(1:len)//'/'//packed_file_name(1:len_1))
      IF (ham_type == 'real') THEN
          Call Input_Packed_Matrices_d(file_directory)
      ELSE IF(ham_type == 'complex') THEN
          Call Input_Packed_Matrices_z(file_directory)
      END IF
      write(iout,*)
      write(iout,*) '          *** Closing Input Matrix File ***'
      CLOSE(50)
      write(iout,*) '          *** Ending Writing of Packed Matrices ***'
      Call IOsys ('close packed_matrices',0,0,0,' ')
      IF( reformat_input_matrix_only) THEN
!            We can opt to quit now if we wish or to re-read in the matrices for
!            a full calculation.
!
          write(iout,*) '          *** Quitting After Reformat ***'
          time(2) = secnds(0.0)
          del_t = time(2) - time(1)
          WRITE(iout,1)
          WRITE(iout,2) del_t
          WRITE(iout,1)
          stop
      ELSE
          write(iout,*)
          write(iout,*) '          *** Reopening File to Read Packed Matrices ***'
          len_1=lenth(packed_file_name)
          Call IOsys('open packed_matrices as old',0,0,0,                            &
                      file_directory(1:len)//'/'//packed_file_name(1:len_1))
          Call Input_Packed_IOsys_Matrix(ham_type,file_directory)
          write(iout,*) '***Closing File Packed Matrices ***'
          Call IOsys('close packed_matrices',0,0,0,' ')
          time(2) = secnds(0.0)
          del_t = time(2) - time(1)
          WRITE(iout,1)
          WRITE(iout,3) del_t
          WRITE(iout,1)
      END IF
  END IF
1 FORMAT('***********************************************'                           &
         '*************************')
2 FORMAT(/,10X,'Time to Reformat/Read/Write/Pack/Cholesky Decompose Matrices = ',f15.8)
3 FORMAT(/,10X,'Time to Open and Read Packed Matrices after Reformat = ',f15.8)
END SUBROUTINE Packed_Matrix
!***********************************************************************
!***********************************************************************
!deck Input_Packed_IOsys_Matrix
!***begin prologue     Input_Packed_IOsys_Matrix
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Hamiltonian input and manipulation.
!***
!***description        This routine can reformat an input Hamiltonian
!                      into a packed form with zeros removed or read and
!                      store already packed matrices for use by the main
!                      routines.  Note that the packed matrices are stored in 
!                      IOsys format for easy reading and writing.
!***references
!***routines called
!***end prologue       Input_Packed_IOsys_Matrix
!
  SUBROUTINE Input_Packed_IOsys_Matrix(ham_type,file_directory)
  IMPLICIT NONE
  CHARACTER(LEN=*)                      :: ham_type
  CHARACTER(LEN=*)                      :: file_directory
  INTEGER                               :: max_buf
  INTEGER                               :: lenth
  INTEGER                               :: len
!
!
  non_zero_overlap_elements = 0
  non_zero_cholesky_elements = 0
  write(iout,*) '          *** Reading in Packed Matrices in Triangle Form ***'
!
!    The matrices here are already packed in column form and only need to be read in
!
  IF (ham_type == 'real') THEN
      IF (non_orth) THEN
          write(iout,*) '          *** Reading Packed Upper Cholesky Factor ***'
          Call IOsys('read integer "number_non_zero_upper_elements" from '//        &
                     'packed_matrices',1,mat_var(1)%number,0,' ')
          ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                &
                   mat_var(1)%row_index(mat_var(1)%number),                         &
                   mat_var(1)%packed_columns_d(mat_var(1)%number),                  &
                   mat_var(1)%matrix_diagonal_d(n3d) )
          Call Read_and_Write_Column_Packed_Matrices                                &
                                           (mat_var(1)%packed_columns_d,            &
                                            mat_var(1)%non_zero_columns,            &
                                            mat_var(1)%row_index,                   &
                                            'read',                                 &
                                            'upper',                                &
                                            mat_var(1)%number,                      &
                                            mat_var(1)%matrix_diagonal_d )
          write(iout,*) '          *** Reading Packed Lower Cholesky Factor ***'
          Call IOsys('read integer "number_non_zero_lower_elements" from '//        &
                     'packed_matrices',1,mat_var(2)%number,0,' ')
          ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                &
                   mat_var(2)%row_index(mat_var(2)%number),                         &
                   mat_var(2)%packed_columns_d(mat_var(2)%number) )
          Call Read_and_Write_Column_Packed_Matrices                                &
                                           (mat_var(2)%packed_columns_d,            &
                                            mat_var(2)%non_zero_columns,            &
                                            mat_var(2)%row_index,                   &
                                            'read',                                 &
                                            'lower',                                &
                                            mat_var(2)%number )
      END IF
      write(iout,*) '          *** Reading Packed Hamiltonian ***'
      Call IOsys('read integer "number_non_zero_hamiltonian_elements" from '//      &
                 'packed_matrices',1,mat_var(3)%number,0,' ')
      ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                    &
               mat_var(3)%row_index(mat_var(3)%number),                             &
               mat_var(3)%packed_columns_d(mat_var(3)%number),                      &
               mat_var(3)%matrix_diagonal_d(n3d) )
      Call Read_and_Write_Column_Packed_Matrices                                    &
                                       (mat_var(3)%packed_columns_d,                &
                                        mat_var(3)%non_zero_columns,                &
                                        mat_var(3)%row_index,                       &
                                       'read',                                      &
                                       'hamiltonian',                               &
                                        mat_var(3)%number,                          &
                                        mat_var(3)%matrix_diagonal_d )
  ELSE IF(ham_type == 'complex') THEN
     IF (non_orth) THEN
         write(iout,*) '          *** Reading Packed Upper Cholesky Factor ***'
         Call IOsys('read integer "number_non_zero_upper_elements" from '//         &
                    'packed_matrices',1,mat_var(1)%number,0,' ')
         ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                 &
                  mat_var(1)%row_index(mat_var(1)%number),                          &
                  mat_var(1)%packed_columns_z(mat_var(1)%number),                   &
                  mat_var(1)%matrix_diagonal_z(n3d) )
         Call Read_and_Write_Column_Packed_Matrices                                 &
                                          (mat_var(1)%packed_columns_z,             &
                                           mat_var(1)%non_zero_columns,             &
                                           mat_var(1)%row_index,                    &
                                           'read',                                  &
                                           'upper',                                 &
                                           mat_var(1)%number,                       &
                                           mat_var(1)%matrix_diagonal_z )
         write(iout,*) '          *** Reading Packed Lower Cholesky Factor ***'
         Call IOsys('read integer "number_non_zero_lower_elements" from '//         &
                    'packed_matrices',1,mat_var(2)%number,0,' ')
         ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                 &
                  mat_var(2)%row_index(mat_var(2)%number),                          &
                  mat_var(2)%packed_columns_z(mat_var(2)%number) )
         Call Read_and_Write_Column_Packed_Matrices                                 &
                                          (mat_var(2)%packed_columns_z,             &
                                           mat_var(2)%non_zero_columns,             &
                                           mat_var(2)%row_index,                    &
                                           'read',                                  &
                                           'lower',                                 &
                                           mat_var(2)%number )
     END IF
     write(iout,*) '          *** Reading Packed Hamiltonian ***'
     Call IOsys('read integer "number_non_zero_hamiltonian_elements" from '//       &
                'packed_matrices',1,mat_var(3)%number,0,' ')
     ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                     &
              mat_var(3)%row_index(mat_var(3)%number),                              &
              mat_var(3)%packed_columns_z(mat_var(3)%number),                       &
              mat_var(3)%matrix_diagonal_z(n3d) )
     Call Read_and_Write_Column_Packed_Matrices                                     &
                                      (mat_var(3)%packed_columns_z,                 &
                                       mat_var(3)%non_zero_columns,                 &
                                       mat_var(3)%row_index,                        &
                                       'read',                                      &
                                       'hamiltonian',                               &
                                       mat_var(3)%number,                           &
                                       mat_var(3)%matrix_diagonal_z )
  END IF
!***********************************************************************
  END SUBROUTINE Input_Packed_IOsys_Matrix
!***********************************************************************
!***********************************************************************
!Deck Input_Packed_Matrices_d
!***begin prologue     Input_Packed_Matrices_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Read in B-spline packed matrices and diagonalize.
!                      Put out information in IOsys format.
!***
!***references
!***routines called
!***end prologue       Input_Packed_Matrices_d
!
  SUBROUTINE Input_Packed_Matrices_d(file_directory)
  IMPLICIT NONE
  REAL*8                                :: gbytes
  INTEGER                               :: len
  INTEGER                               :: len_1
  INTEGER                               :: lenth
  CHARACTER(LEN=*)                      :: file_directory
  CHARACTER(LEN=80)                     :: matrix_file_name
!
  gbytes = tri_size*8*3
!
!                     Did not comment this routine extensively since its
!                     Hermitian version of the previous routine.
!
! This option is just used for test purposes and if this option is set, we exit 
! immediately after the diagonalization.
!
  len = lenth(file_directory)
  IF (diagonalize_only) THEN
      Call Diagonalize_Packed_Matrices_d('disk',                                       &
                                         n3d,                                          &
                                         .true.)
      stop
  END IF
!
! This is the option to reformat the input matrix in IOsys format and then return 
! to the calling routine.

!
!    Triangular Packing by Columns.
!
  IF(non_orth) THEN      
     ALLOCATE(triangle_overlap_d(tri_size))
     Call Fill_Matrix_from_Disk (                                                       &
                                 triangle_overlap_d,                                    &
                                 non_zero_overlap_elements,                             &
                                 50)
     write(iout,*) 'Number Non_Zero Overlap In = ',non_zero_overlap_elements
     IF (overlap_to_disk) THEN   
         drop_tol=drop_overlap
         write(iout,*)
         write(iout,*) '          *** Pack the Overlap ***'
         ALLOCATE(mat_var(4)%non_zero_columns(n3d),                                     &
                  mat_var(4)%row_index(lenbuf),                                         &
                  mat_var(4)%packed_columns_d(lenbuf) ,                                 &
                  mat_var(4)%matrix_diagonal_d(n3d) )
         Call Pack_Triangle ('overlap',                                                 &
                              triangle_overlap_d,                                       &
                              mat_var(4)%packed_columns_d,                              &
                              mat_var(4)%non_zero_columns,                              &
                              mat_var(4)%row_index,                                     &
                              mat_var(4)%number,                                        &
                              mat_var(4)%matrix_diagonal_d)
         DEALLOCATE(mat_var(4)%non_zero_columns,                                        &
                    mat_var(4)%row_index,                                               &
                    mat_var(4)%packed_columns_d,                                        &
                    mat_var(4)%matrix_diagonal_d )
     END IF
     ALLOCATE(upper_d(tri_size))
     upper_d(1:tri_size) = triangle_overlap_d(1:tri_size)
     write(iout,*)
     write(iout,*)'          ***Form Cholesky Decomposition ***'
     call dpptrf('u',n3d,upper_d,info)
     IF(info /= 0) THEN
        write(iout,*) '          *** Singular Overlap Matrix:Quit ***'
        Call lnkerr('Quit.  Singular Overlap Matrix')
     END IF
     write(iout,*)
     write(iout,*) '          *** Pack the Cholesky Factors ***'
     drop_tol=drop_cholesky
     ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                         &
              mat_var(1)%row_index(lenbuf),                                             &
              mat_var(1)%packed_columns_d(lenbuf) ,                                     &
              mat_var(1)%matrix_diagonal_d(n3d) )
     Call Pack_Triangle ('upper',                                                       &
                          upper_d,                                                      &
                          mat_var(1)%packed_columns_d,                                  &
                          mat_var(1)%non_zero_columns,                                  &
                          mat_var(1)%row_index,                                         &
                          mat_var(1)%number,                                            &
                          mat_var(1)%matrix_diagonal_d)
     DEALLOCATE(mat_var(1)%non_zero_columns,                                            &
                mat_var(1)%row_index,                                                   &
                mat_var(1)%packed_columns_d )
     Call Upper_to_Lower(upper_d,triangle_overlap_d)
     DEALLOCATE(upper_d)
     ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                         &
              mat_var(2)%row_index(lenbuf),                                             &
              mat_var(2)%packed_columns_d(lenbuf) )
     Call Pack_Triangle ('lower',                                                       &
                          triangle_overlap_d,                                           &
                          mat_var(2)%packed_columns_d,                                  &
                          mat_var(2)%non_zero_columns,                                  &
                          mat_var(2)%row_index,                                         &
                          mat_var(2)%number,                                            &
                          mat_var(1)%matrix_diagonal_d)
     DEALLOCATE(mat_var(2)%non_zero_columns,                                            &
                mat_var(2)%row_index,                                                   &
                mat_var(2)%packed_columns_d,                                            &
                mat_var(1)%matrix_diagonal_d )
     DEALLOCATE(triangle_overlap_d )
  END IF
  ALLOCATE(triangle_hamiltonian_d(1:tri_size))
  Call Fill_Matrix_from_Disk (                                                          &
                              triangle_hamiltonian_d,                                   &
                              non_zero_hamiltonian_elements,                            &
                              50)
  write(iout,*) 'Number Non_Zero Hamiltonian In = ',non_zero_hamiltonian_elements
  write(iout,*)
  write(iout,*) '          *** Pack the Hamiltonian ***'
  ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                            &
           mat_var(3)%row_index(lenbuf),                                                &
           mat_var(3)%packed_columns_d(lenbuf),                                         &
           mat_var(3)%matrix_diagonal_d(n3d) )
  drop_tol=drop_hamiltonian
  Call Pack_Triangle ('hamiltonian',                                                    &
                       triangle_hamiltonian_d,                                          &
                       mat_var(3)%packed_columns_d,                                     &
                       mat_var(3)%non_zero_columns,                                     &
                       mat_var(3)%row_index,                                            &
                       mat_var(3)%number,                                               &
                       mat_var(3)%matrix_diagonal_d)
  DEALLOCATE(mat_var(3)%non_zero_columns,                                               &
             mat_var(3)%row_index,                                                      &
             mat_var(3)%packed_columns_d,                                               &
             mat_var(3)%matrix_diagonal_d )
  DEALLOCATE(triangle_hamiltonian_d)
!***********************************************************************
  END SUBROUTINE Input_Packed_Matrices_d
!***********************************************************************
!***********************************************************************
!Deck Input_Packed_Matrices_z
!***begin prologue     Input_Packed_Matrices_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Read in B-spline packed matrices and diagonalize.
!                      Put out information in IOsys format.
!***
!***references
!***routines called
!***end prologue       Input_Packed_Matrices_z
!
  SUBROUTINE Input_Packed_Matrices_z(file_directory)
  IMPLICIT NONE
  REAL*8                                :: gbytes
  INTEGER                               :: len
  INTEGER                               :: len_1
  INTEGER                               :: lenth
  CHARACTER(LEN=*)                      :: file_directory
  CHARACTER(LEN=80)                     :: matrix_file_name
  gbytes = tri_size*16*3
!
!                     Did not comment this routine extensively since its
!                     Hermitian version of the previous routine.
!
! This option is just used for test purposes and if this option is set, we exit 
! immediately after the diagonalization.
!
  len = lenth(file_directory)
  IF (diagonalize_only) THEN
      Call Diagonalize_Packed_Matrices_z('disk',                                            &
                                         n3d,                                               &
                                         .false.)
      stop
  END IF
!
! This is the option to reformat the input matrix in IOsys format and then return 
! to the calling routine.

!
!       Triangular Packing by Columns.
!
  IF(non_orth) THEN      
     ALLOCATE(triangle_overlap_z(tri_size))
     Call Fill_Matrix_from_Disk (                                                        &
                                 triangle_overlap_z,                                     &
                                 non_zero_overlap_elements,                              &
                                 50)
     write(iout,*) 'Number Non_Zero Overlap In = ',non_zero_overlap_elements
     IF (overlap_to_disk) THEN   
         drop_tol=drop_overlap
         write(iout,*)
         write(iout,*) 'Pack the Overlap'
         ALLOCATE(mat_var(4)%non_zero_columns(n3d),                                      &
                  mat_var(4)%row_index(lenbuf),                                          &
                  mat_var(4)%packed_columns_z(lenbuf) ,                                  &
                  mat_var(4)%matrix_diagonal_z(n3d) )
         Call Pack_Triangle ('overlap',                                                  &
                              triangle_overlap_z,                                        &
                              mat_var(4)%packed_columns_z,                               &
                              mat_var(4)%non_zero_columns,                               &
                              mat_var(4)%row_index,                                      &
                              mat_var(4)%number,                                         &
                              mat_var(4)%matrix_diagonal_z)    
         DEALLOCATE(mat_var(4)%non_zero_columns,                                         &
                    mat_var(4)%row_index,                                                &
                    mat_var(4)%packed_columns_z,                                         &
                    mat_var(4)%matrix_diagonal_z )
     END IF
     ALLOCATE(upper_z(tri_size))
     upper_z(1:tri_size) = triangle_overlap_z(1:tri_size)
     write(iout,*)
     write(iout,*)'          *** Form Cholesky Decomposition ***'
     call dpptrf('u',n3d,upper_z,info)
     IF(info /= 0) THEN
        write(iout,*) '          *** Singular Overlap Matrix:Quit ***'
        Call lnkerr('Quit.  Singular Overlap Matrix')
     END IF
     write(iout,*)
     write(iout,*) '          *** Pack the Cholesky Factors ***'
     drop_tol=drop_cholesky
     ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                          &
              mat_var(1)%row_index(lenbuf),                                              &
              mat_var(1)%packed_columns_z(lenbuf) ,                                      &
              mat_var(1)%matrix_diagonal_z(n3d) )
     Call Pack_Triangle ('upper',                                                        &
                          upper_z,                                                       &
                          mat_var(1)%packed_columns_z,                                   &
                          mat_var(1)%non_zero_columns,                                   &
                          mat_var(1)%row_index,                                          &
                          mat_var(1)%number,                                             &
                          mat_var(1)%matrix_diagonal_z)
     DEALLOCATE(mat_var(1)%non_zero_columns,                                             &
                mat_var(1)%row_index,                                                    &
                mat_var(1)%packed_columns_z )
     Call Upper_to_Lower(upper_z,triangle_overlap_z)
     DEALLOCATE(upper_z)
     ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                          &
              mat_var(2)%row_index(lenbuf),                                              &
              mat_var(2)%packed_columns_z(lenbuf) )
     Call Pack_Triangle ('lower',                                                        &
                          triangle_overlap_z,                                            &
                          mat_var(2)%packed_columns_z,                                   &
                          mat_var(2)%non_zero_columns,                                   &
                          mat_var(2)%row_index,                                          &
                          mat_var(2)%number,                                             &
                          mat_var(1)%matrix_diagonal_z)
     DEALLOCATE(mat_var(2)%non_zero_columns,                                             &
                mat_var(2)%row_index,                                                    &
                mat_var(2)%packed_columns_z,                                             &
                mat_var(1)%matrix_diagonal_z )
     DEALLOCATE(triangle_overlap_z )
  END IF
  ALLOCATE(triangle_hamiltonian_z(1:tri_size))
  Call Fill_Matrix_from_Disk (                                                           &
                              triangle_hamiltonian_z,                                    &
                              non_zero_hamiltonian_elements,                             &
                              50)
  write(iout,*) 'Number Non_Zero Hamiltonian In = ',non_zero_hamiltonian_elements
  write(iout,*)
  write(iout,*) '          *** Pack the Hamiltonian ***'
  ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                             &
           mat_var(3)%row_index(lenbuf),                                                 &
           mat_var(3)%packed_columns_z(lenbuf),                                          &
           mat_var(3)%matrix_diagonal_z(n3d) )
  drop_tol=drop_hamiltonian
  Call Pack_Triangle ('hamiltonian',                                                     &
                       triangle_overlap_z,                                               &
                       mat_var(3)%packed_columns_z,                                      &
                       mat_var(3)%non_zero_columns,                                      &
                       mat_var(3)%row_index,                                             &
                       mat_var(3)%number,                                                &
                       mat_var(3)%matrix_diagonal_z)
  DEALLOCATE(mat_var(3)%non_zero_columns,                                                &
             mat_var(3)%row_index,                                                       &
             mat_var(3)%packed_columns_z,                                                &
             mat_var(3)%matrix_diagonal_z )
  DEALLOCATE(triangle_hamiltonian_z)
!***********************************************************************
  END SUBROUTINE Input_Packed_Matrices_z
!***********************************************************************
!***********************************************************************
!Deck Diagonalize_Packed_Matrices_d
!***begin prologue     Diagonalize_Packed_Matrices_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Read in B-spline packed matrices and diagonalize.
!                      Put out information in IOsys format.
!***
!***references
!***routines called
!***end prologue       Diagonalize_Packed_Matrices_d
!
  SUBROUTINE Diagonalize_Packed_Matrices_d (matrix_source,mat_size,                          &
                                            get_eigenvectors,to_disk,                        &
                                            file_key,ham,s)
  IMPLICIT NONE
  CHARACTER(LEN=*)                          :: matrix_source
  LOGICAL                                   :: get_eigenvectors 
  INTEGER                                   :: mat_size
  LOGICAL, OPTIONAL                         :: to_disk 
  CHARACTER (LEN=*), OPTIONAL               :: file_key 
  REAL*8, DIMENSION(:), OPTIONAL            :: ham
  REAL*8, DIMENSION(:), OPTIONAL            :: s
  CHARACTER(LEN=1)                          :: n_v
  CHARACTER*80                              :: overlap 
  CHARACTER*80                              :: hamiltonian 
  CHARACTER*8                               :: keywrd
  IF ( present(file_key) == .false. ) THEN
       overlap = 'overlap'
       hamiltonian = 'hamiltonian'
  ELSE
       overlap = 'overlap_matrix_'//file_key
       hamiltonian = 'hamiltonian_matrix_'//file_key
  END IF
  lwork = 3 * mat_size
  IF (non_orth) THEN
      lwork = 10 * mat_size
  END IF
  ALLOCATE( eig(1:mat_size), work_d(1:lwork) )
  IF ( matrix_source /= 'input') THEN
       ALLOCATE( triangle_hamiltonian_d(1:tri_size) )
       IF (non_orth) THEN
           ALLOCATE( triangle_overlap_d(1:tri_size) )
       END IF
  END IF
  n_v='n'
  IF (get_eigenvectors == .true. ) THEN
      ALLOCATE(eigenvectors_d(1:mat_size,1:mat_size))
      n_v='v'
  END IF
  IF (matrix_source == 'disk' ) THEN
      write(iout,*) '          *** Diagonalizing Big Matrix From Disk ***'
      IF(non_orth) THEN
         Call Fill_Matrix_from_Disk(triangle_overlap_d,                                     &
                                    non_zero_overlap_elements,                              &
                                    50)
         Call Fill_Matrix_from_Disk(triangle_hamiltonian_d,                                 &
                                    non_zero_hamiltonian_elements,                          &
                                    50)
         call dspgv(1,n_v,'u',mat_size,triangle_hamiltonian_d,triangle_overlap_d,           &
                                  eig,eigenvectors_d,mat_size,work_d,lwork,info)
      ELSE 
         Call Fill_Matrix_from_Disk(triangle_hamiltonian_d,                                 &
                                    non_zero_hamiltonian_elements,                          &
                                    50)
         call dspev(n_v,'u',mat_size,triangle_hamiltonian_d,eig,eigenvectors_d,mat_size,work_d,info)
      END IF
  ELSE IF(matrix_source == 'buffers' ) THEN
      write(iout,*) '          *** Diagonalizing Big Matrix From Buffers ***'
      IF (non_orth) THEN
          Call IOsys('read integer number_non_zero_overlap_elements from '//                &
                     'packed_matrices',1,mat_var(1)%number,0,' ')
          ALLOCATE(mat_var(1)%non_zero_columns(mat_size),                                   &
                   mat_var(1)%row_index(mat_var(1)%number),                                 &
                   mat_var(1)%packed_columns_d(mat_var(1)%number),                          &
                   mat_var(1)%matrix_diagonal_d(mat_size) )
          Call Read_and_Write_Column_Packed_Matrices (                                      &
                                          mat_var(1)%packed_columns_d,                      &
                                          mat_var(1)%non_zero_columns,                      &
                                          mat_var(1)%row_index,                             &
                                          'read',                                           &
                                          overlap,                                          &
                                          mat_var(1)%number,                                &
                                          mat_var(1)%matrix_diagonal_d )
          Call Fill_Upper_Triangular_Matrix_from_Column_Buffers_d (                         &
                                          triangle_overlap_d,                               &
                                          mat_var(1)%matrix_diagonal_d,                     &
                                          mat_var(1)%packed_columns_d,                      &
                                          mat_var(1)%non_zero_columns,                      &
                                          mat_var(1)%row_index)
          DEALLOCATE(mat_var(1)%row_index, mat_var(1)%packed_columns_d )
          Call IOsys('read integer number_non_zero_hamiltonian_elements from '//            &
                     'packed_matrices',1,mat_var(1)%number,0,' ') 
          ALLOCATE(mat_var(1)%row_index(mat_var(1)%number),                                 &
                   mat_var(1)%packed_columns_d(mat_var(1)%number) )
          Call Read_and_Write_Column_Packed_Matrices (                                      &
                                            mat_var(1)%packed_columns_d,                    &
                                            mat_var(1)%non_zero_columns,                    &
                                            mat_var(1)%row_index,                           &
                                            'read',                                         &
                                            hamiltonian,                                    &
                                            mat_var(1)%number,                              &
                                            mat_var(1)%matrix_diagonal_d )
          Call Fill_Upper_Triangular_Matrix_from_Column_Buffers_d (                         &
                                              triangle_hamiltonian_d,                       &
                                              mat_var(1)%matrix_diagonal_d,                 &
                                              mat_var(1)%packed_columns_d,                  &
                                              mat_var(1)%non_zero_columns,                  &
                                              mat_var(1)%row_index)
          call dspgv(1,n_v,'u',mat_size,triangle_hamiltonian_d,                             &
                                  triangle_overlap_d,                                       &
                                  eig,eigenvectors_d,mat_size,work_d,lwork,info)
      ELSE
          Call IOsys('read integer number_non_zero_hamiltonian_elements from '//            &
                     'packed_matrices',1,mat_var(1)%number,0,' ') 
          ALLOCATE(mat_var(1)%non_zero_columns(mat_size),                                   &
                   mat_var(1)%row_index(mat_var(1)%number),                                 &
                   mat_var(1)%packed_columns_d(mat_var(1)%number),                          &
                   mat_var(1)%matrix_diagonal_d(mat_size) )
          Call Read_and_Write_Column_Packed_Matrices (                                      &
                                            mat_var(1)%packed_columns_d,                    &
                                            mat_var(1)%non_zero_columns,                    &
                                            mat_var(1)%row_index,                           &
                                            'read',                                         &
                                            hamiltonian,                                    &
                                            mat_var(1)%number,                              &
                                            mat_var(1)%matrix_diagonal_d )
          Call Fill_Upper_Triangular_Matrix_from_Column_Buffers_d (                         &
                                              triangle_hamiltonian_d,                       &
                                              mat_var(1)%matrix_diagonal_d,                 &
                                              mat_var(1)%packed_columns_d,                  &
                                              mat_var(1)%non_zero_columns,                  &
                                              mat_var(1)%row_index)
          call dspev(n_v,'u',mat_size,triangle_hamiltonian_d,eig,eigenvectors_d,            &
                     mat_size,work_d,info)
      END IF
      DEALLOCATE(mat_var(1)%non_zero_columns,                                               &
                 mat_var(1)%row_index,                                                      &
                 mat_var(1)%packed_columns_d,                                               &
                 mat_var(1)%matrix_diagonal_d )
  ELSE IF(matrix_source == 'input') THEN
      write(iout,*) '          *** Diagonalizing Big Matrix From Input Matrices ***'
      IF (non_orth) THEN
          call dspgv(1,n_v,'u',mat_size,ham,s,eig,eigenvectors_d,mat_size,work_d,lwork,info)
      ELSE
          call dspev(n_v,'u',mat_size,ham,eig,eigenvectors_d,mat_size,work_d,info)
      END IF 
  END IF
  Call IOsys('write real "hamiltonian_eigenvalues_'//file_key//'" to '//'packed_matrices',  &
              mat_size,eig,0,' ')
  IF (get_eigenvectors == .true. ) THEN
      Call IOsys('write real "hamiltonian_eigenvectors_'//file_key//'" to '//'packed_matrices',  &
                  mat_size*mat_size,eigenvectors_d,0,' ')
      write(iout,*) 'Writing Eigenvector'
      keywrd='real'
      REWIND(20)
      WRITE(20) keywrd
      WRITE(20) eigenvectors_d(:,1)
      CLOSE(20)
      DEALLOCATE(eigenvectors_d)
  END IF
  local_title='eigenvalues'
  call prntfmn(local_title,eig,eigenvectors_to_print,1,eigenvectors_to_print,1,iout,'e')
  write(iout,*) '          *** Done Diagonalizing Big Matrix for this Symmetry ***'
  DEALLOCATE( eig, work_d )
  IF ( matrix_source /= 'input') THEN
       DEALLOCATE( triangle_hamiltonian_d )
       IF (non_orth) THEN
           DEALLOCATE( triangle_overlap_d )
       END IF
  END IF
!***********************************************************************
  END SUBROUTINE Diagonalize_Packed_Matrices_d
!***********************************************************************
!***********************************************************************
!Deck Diagonalize_Packed_Matrices_z
!***begin prologue     Diagonalize_Packed_Matrices_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Read in B-spline packed matrices and diagonalize.
!                      Put out information in IOsys format.
!***
!***references
!***routines called
!***end prologue       Diagonalize_Packed_Matrices_z
!
  SUBROUTINE Diagonalize_Packed_Matrices_z (matrix_source,mat_size,get_eigenvectors,     &
                                            to_disk,file_key,ham,s)
  IMPLICIT NONE
  CHARACTER(LEN=*)                          :: matrix_source
  LOGICAL                                   :: get_eigenvectors 
  INTEGER                                   :: mat_size
  LOGICAL, OPTIONAL                         :: to_disk 
  CHARACTER(LEN=*), OPTIONAL                :: file_key 
  COMPLEX*16, DIMENSION(:), OPTIONAL        :: ham
  COMPLEX*16, DIMENSION(:), OPTIONAL        :: s
  CHARACTER(LEN=1)                          :: n_v
  CHARACTER(LEN=80)                         :: overlap
  CHARACTER(LEN=80)                         :: hamiltonian
  CHARACTER*8                               :: keywrd
  IF ( present(file_key) == .false. ) THEN
       overlap = 'overlap'
       hamiltonian = 'hamiltonian'
  ELSE
       overlap = 'overlap_matrix_'//file_key
       hamiltonian = 'hamiltonian_matrix_'//file_key
  END IF
  lwork = 3 * mat_size
  IF (non_orth) THEN
      lwork = 10 * mat_size
  END IF
  ALLOCATE( eig(1:mat_size), work_z(1:lwork), rwork(1:lwork) )
  IF ( matrix_source /= 'input') THEN
       ALLOCATE( triangle_hamiltonian_z(1:tri_size) )
       IF (non_orth) THEN
           ALLOCATE( triangle_overlap_z(1:tri_size) )
       END IF
  END IF
  n_v='n'
  IF (get_eigenvectors == .true. ) THEN
      ALLOCATE(eigenvectors_z(1:mat_size,1:mat_size))
      n_v='v'
  END IF
  IF (matrix_source == 'disk' ) THEN
      write(iout,*) '          *** Diagonalizing Big Matrix From Disk ***'
      IF(non_orth) THEN
         Call Fill_Matrix_from_Disk(triangle_overlap_z,                                     &
                                    non_zero_overlap_elements,                              &
                                    50)
         Call Fill_Matrix_from_Disk(triangle_hamiltonian_z,                                 &
                                    non_zero_hamiltonian_elements,                          &
                                    50)
         call zhpgv(1,n_v,'u',mat_size,triangle_hamiltonian_z,                              &
                                          triangle_overlap_z,                               &
                                          eig,eigenvectors_z,mat_size,work_z,rwork,info)
      ELSE 
         Call Fill_Matrix_from_Disk(triangle_hamiltonian_z,                                 &
                                    non_zero_hamiltonian_elements,                          &
                                    50)
         call zhpev(n_v,'u',mat_size,triangle_hamiltonian_z,eig,eigenvectors_z,             &
                            mat_size,work_z,rwork,info)
      END IF
  ELSE IF(matrix_source == 'buffers' ) THEN
      write(iout,*) '          *** Diagonalizing Big Matrix From Buffers ***'
      Call IOsys('read integer number_non_zero_overlap_elements from '//                    &
                 'packed_matrices',1,mat_var(1)%number,0,' ')
      ALLOCATE(mat_var(1)%non_zero_columns(mat_size),                                   &
               mat_var(1)%row_index(mat_var(1)%number),                                 &
               mat_var(1)%packed_columns_z(mat_var(1)%number),                          &
               mat_var(1)%matrix_diagonal_z(mat_size) )
      IF (non_orth) THEN
          Call Read_and_Write_Column_Packed_Matrices (                                      &
                                          mat_var(1)%packed_columns_z,                      &
                                          mat_var(1)%non_zero_columns,                      &
                                          mat_var(1)%row_index,                             &
                                          'read',                                           &
                                          overlap,                                          &
                                          mat_var(1)%number,                                &
                                          mat_var(1)%matrix_diagonal_z )
          Call Fill_Upper_Triangular_Matrix_from_Column_Buffers_z (                         &
                                          triangle_overlap_z,                               &
                                          mat_var(1)%matrix_diagonal_z,                     &
                                          mat_var(1)%packed_columns_z,                      &
                                          mat_var(1)%non_zero_columns,                      &
                                          mat_var(1)%row_index)
          DEALLOCATE(mat_var(1)%row_index, mat_var(1)%packed_columns_z )
          Call IOsys('read integer number_non_zero_hamiltonian_elements from '//            &
                     'packed_matrices',1,mat_var(1)%number,0,' ') 
          ALLOCATE(mat_var(1)%row_index(mat_var(1)%number),                                 &
                   mat_var(1)%packed_columns_z(mat_var(1)%number) )
          Call Read_and_Write_Column_Packed_Matrices (                                      &
                                            mat_var(1)%packed_columns_z,                    &
                                            mat_var(1)%non_zero_columns,                    &
                                            mat_var(1)%row_index,                           &
                                            'read',                                         &
                                            hamiltonian,                                    &
                                            mat_var(1)%number,                              &
                                            mat_var(1)%matrix_diagonal_z )
          Call Fill_Upper_Triangular_Matrix_from_Column_Buffers_z (                         &
                                              triangle_hamiltonian_z,                       &
                                              mat_var(1)%matrix_diagonal_z,                 &
                                              mat_var(1)%packed_columns_z,                  &
                                              mat_var(1)%non_zero_columns,                  &
                                              mat_var(1)%row_index)
          call zhpgv(1,n_v,'u',mat_size,triangle_hamiltonian_z,                             &
                                           triangle_overlap_z,                              &
                                           eig,eigenvectors_z,mat_size,work_z,              &
                                           rwork,info)
      ELSE
          Call IOsys('read integer number_non_zero_hamiltonian_elements from '//            &
                     'packed_matrices',1,mat_var(1)%number,0,' ') 
          ALLOCATE(mat_var(1)%non_zero_columns(mat_size),                                   &
                   mat_var(1)%row_index(mat_var(1)%number),                                 &
                   mat_var(1)%packed_columns_z(mat_var(1)%number),                          &
                   mat_var(1)%matrix_diagonal_z(mat_size) )
          Call Read_and_Write_Column_Packed_Matrices (                                      &
                                            mat_var(1)%packed_columns_z,                    &
                                            mat_var(1)%non_zero_columns,                    &
                                            mat_var(1)%row_index,                           &
                                            'read',                                         &
                                            hamiltonian,                                    &
                                            mat_var(1)%number,                              &
                                            mat_var(1)%matrix_diagonal_z )
          Call Fill_Upper_Triangular_Matrix_from_Column_Buffers_z (                         &
                                              triangle_hamiltonian_z,                       &
                                              mat_var(1)%matrix_diagonal_z,                 &
                                              mat_var(1)%packed_columns_z,                  &
                                              mat_var(1)%non_zero_columns,                  &
                                              mat_var(1)%row_index)
          call zhpev(n_v,'u',mat_size,triangle_hamiltonian_z,                               &
                             eig,eigenvectors_z,mat_size,work_z,rwork,info)
      END IF
      DEALLOCATE(mat_var(1)%non_zero_columns,                                               &
                 mat_var(1)%row_index,                                                      &
                 mat_var(1)%packed_columns_z,                                               &
                 mat_var(1)%matrix_diagonal_z )
  ELSE IF(matrix_source == 'input') THEN
      write(iout,*) '          *** Diagonalizing Big Matrix From Input Matrices ***'
      IF (non_orth) THEN
         call zhpgv(1,n_v,'u',mat_size,ham,                                                 &
                                          s,                                                &
                                          eig,                                              &
                                          eigenvectors_z,                                   &
                                          mat_size,                                         &
                                          work_z,                                           &
                                          rwork,                                            &
                                          info)
      ELSE
          call zhpev(n_v,'u',mat_size,ham,                                                  &
                                         eig,                                               &
                                         eigenvectors_z,                                    &
                                         mat_size,                                          &
                                         work_z,                                            &
                                         rwork,                                             &
                                         info)
      END IF 
  END IF
  Call IOsys('write real "hamiltonian_eigenvalues_'//file_key//'" to '//'packed_matrices',  &
              mat_size,eig,0,' ')
  IF (get_eigenvectors == .true. ) THEN
      Call IOsys('write real "hamiltonian_eigenvectors_'//file_key//'" to '//'packed_matrices',  &
                  2*mat_size*mat_size,eigenvectors_z,0,' ')
      keywrd='complex'
      REWIND(20)
      WRITE(20) keywrd
      WRITE(20) eigenvectors_z(:,1)
      CLOSE(20)
      DEALLOCATE(eigenvectors_z)
  END IF
  local_title='eigenvalues'
  call prntfmn(local_title,eig,eigenvectors_to_print,1,eigenvectors_to_print,1,iout,'e')
  write(iout,*) '          *** Done Diagonalizing Big Matrix for this Symmetry ***'
  DEALLOCATE( eig, work_z )
  IF ( matrix_source /= 'input') THEN
       DEALLOCATE( triangle_hamiltonian_z )
       IF (non_orth) THEN
           DEALLOCATE( triangle_overlap_z )
       END IF
  END IF
!***********************************************************************
  END SUBROUTINE Diagonalize_Packed_Matrices_z
!***********************************************************************
!***********************************************************************
  END  MODULE Packed_Matrix_Module
!***********************************************************************
!***********************************************************************
