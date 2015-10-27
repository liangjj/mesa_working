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
  write(iout,*)
  IF (input_matrices == 'packed_in_drake_format') THEN
      write(iout,*) '          Reformatting Packed Matrices to IOsys Format'
!
!     In this branch of the if, we come in with a packed matrix from the Drake codes, 
!     and go out with a packed matrix in IOsys format.  The packed matrix from IOsys 
!     is done in the standard or triangular manner.  In standard manner an array of 
!     two indices and an array of elements that correspond to those 
!     indices is output along with the value of the number of non zero elements.  
!     Only the lower triangle is packed.  The other alternative is a packing of the 
!     triangle of the matrix.  This column oriented packing requires more extensive
!     logic but is needed to efficiently solve the linear equations when a 
!     non-orthogonal basis is used.
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
          Call Input_Packed_Matrices_d(file_directory,output_matrices)
      ELSE IF(ham_type == 'complex') THEN
          Call Input_Packed_Matrices_z(file_directory,output_matrices)
      END IF
      write(iout,*)
      write(iout,*) '          *** Closing Input Matrix File ***'
      CLOSE(50)
      write(iout,*) '          *** Ending Writing of Packed Matrices ***'
      Call IOsys ('close packed_matrices',0,0,0,' ')
      IF( reformat_input_matrix_only) THEN
!
!         We can opt to quit now if we wish or to re-read in the matrices for
!         a full calculation.
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
          Call Input_Packed_IOsys_Matrix(ham_type,file_directory,output_matrices)
          write(iout,*) '***Closing File Packed Matrices ***'
          Call IOsys('close packed_matrices',0,0,0,' ')
          time(2) = secnds(0.0)
          del_t = time(2) - time(1)
          WRITE(iout,1)
          WRITE(iout,3) del_t
          WRITE(iout,1)
      END IF
  ELSE IF(input_matrices =='packed_in_iosys_format') THEN
!
!         Take the already IOsys formatted files and read in the matrices.
!
          time(1) = secnds(0.0)
          write(iout,*) '          Reading in the IOsys Formated Matrices'
          len=lenth(file_directory)
          write(iout,*) '          *** Opening File to Read Packed Matrices ***'
          len_1=lenth(packed_file_name)
          Call IOsys('open packed_matrices as old',0,0,0,                            &
                      file_directory(1:len)//'/'//packed_file_name(1:len_1))
          Call Input_Packed_IOsys_Matrix(ham_type,file_directory,output_matrices)
          write(iout,*) '          *** Closing File Packed Matrices ***'
          Call IOsys('close packed_matrices',0,0,0,' ')
          time(2) = secnds(0.0)
          del_t = time(2) - time(1)
          WRITE(iout,1)
          WRITE(iout,4) del_t
          WRITE(iout,1)
  ELSE
          Write(iout,5) input_matrices
          Call lnkerr('Bad Input Matrix Keyword')
  END IF
1 FORMAT('***********************************************'                           &
         '*************************')
2 FORMAT(/,10X,'Time to Reformat/Read/Write/Pack/Cholesky Decompose Matrices = ',f15.8)
3 FORMAT(/,10X,'Time to Open and Read Packed Matrices after Reformat = ',f15.8)
4 FORMAT(/,10X,'Time to Open and Read Packed Matrices = ',f15.8)
5 FORMAT(/,10X,'Bad Input Matrix Keyword Used = ',a80,/,10x,                         &
               'Allowed Values are: packed_in_drake_format',/,10x,                   &
               '                    packed_in_iosys_format')
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
  SUBROUTINE Input_Packed_IOsys_Matrix(ham_type,file_directory,output_matrices)
  IMPLICIT NONE
  CHARACTER(LEN=*)                      :: ham_type
  CHARACTER(LEN=*)                      :: file_directory
  CHARACTER(LEN=*)                      :: output_matrices
  INTEGER                               :: max_buf
  INTEGER                               :: lenth
  INTEGER                               :: len
!
!
  non_zero_overlap_elements = 0
  non_zero_cholesky_elements = 0
  IF(output_matrices == 'packed_in_standard_iosys_format') THEN
!
!    If the matrices were written using a two buffer packed representation
!    this section of code will read them in and reformat them into the packed
!    column representation needed in the code.  With the exception of the overlap
!    matrix, which is no longer needed once the Cholesky factors are available,
!    the required arrays remain allocated after being filled and are available.
!
     write(iout,*) '          ***Reading in Packed Matrices in Standard Two Buffer'//  &
                   ' Form ***'
     Call IOsys('read integer number_non_zero_hamiltonian_elements from '//            &
                'packed_matrices',1,non_zero_hamiltonian_elements,0,' ') 
     IF (non_orth) THEN
         Call IOsys('read integer number_non_zero_overlap_elements from '//            &
                    'packed_matrices',1,non_zero_overlap_elements,0,' ')
     END IF
     max_buf = max(non_zero_overlap_elements,non_zero_hamiltonian_elements)
     ALLOCATE(ibuf(2,max_buf), h_buf_d(max_buf), diag_d(n3d))
!
!          In this section we read in the already packed matrices and then
!          fill up the full matrices with the elements.
!          
     IF (ham_type == 'real') THEN
         IF (non_orth) THEN
!
             IF (overlap_to_disk) THEN
                 ALLOCATE(triangle_overlap_d(tri_size))
!
!                Read in the packed matrices into the buffers.
!
                 Call Read_and_Write_Packed_Matrices(ibuf,                             &
                                                     h_buf_d,                          &
                                                    'read',                            &
                                                    'overlap',                         &
                                                     non_zero_overlap_elements,        &
                                                     diag_d )
!
!                Fill the full matrix from the buffers.
!
                 Call Fill_Matrix_from_Buffers(triangle_overlap_d,                     &
                                               ibuf,                                   &
                                               h_buf_d,                                &
                                               diag_d,                                 &
                                               non_zero_overlap_elements)
                 write(iout,*) '          *** Overlap Read in.  Non-Zero Elements *** '&
                                                   '= ', non_zero_overlap_elements
!
!                Form the column packed representation and write it out to disk.
!
                 drop_tol=drop_overlap
                 write(iout,*)
                 write(iout,*) '          *** Pack the Overlap ***'
                 ALLOCATE(mat_var(4)%non_zero_columns(n3d),                            &
                          mat_var(4)%row_index(lenbuf),                                &
                          mat_var(4)%packed_columns_d(lenbuf) ,                        &
                          mat_var(4)%matrix_diagonal_d(n3d) )
                 Call Pack_Triangle ('overlap',                                        &
                                     triangle_overlap_d,                               &
                                     mat_var(4)%packed_columns_d,                      &
                                     mat_var(4)%non_zero_columns,                      &
                                     mat_var(4)%row_index,                             &
                                     mat_var(4)%number,                                &
                                     mat_var(4)%matrix_diagonal_d)
                 DEALLOCATE(mat_var(4)%non_zero_columns,                               &
                            mat_var(4)%row_index,                                      &
                            mat_var(4)%packed_columns_d,                               &
                            mat_var(4)%matrix_diagonal_d )
                 DEALLOCATE(triangle_overlap_d)
             END IF
!
!
!            The Cholesky factors have already been packed.  Read them in.
!
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
!
!        Do similarly for Hamiltonian
!
         ALLOCATE(triangle_hamiltonian_d(tri_size))
         Call Read_and_Write_Packed_Matrices(ibuf,                                     &
                                             h_buf_d,                                  &
                                             'read',                                   &
                                             'hamiltonian',                            &
                                             non_zero_hamiltonian_elements,            & 
                                             diag_d )
         Call Fill_Matrix_from_Buffers(triangle_hamiltonian_d,                         &
                                       ibuf,                                           &
                                       h_buf_d,                                        &
                                       diag_d,                                         &
                                       non_zero_hamiltonian_elements)
         write(iout,*) '          *** Hamiltonian Read in.  Non-Zero Elements *** = ', &
                                                            non_zero_hamiltonian_elements
         write(iout,*)
         write(iout,*) '          *** Pack the Hamiltonian ***'
         ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                    &
                  mat_var(3)%row_index(lenbuf),                                        &
                  mat_var(3)%packed_columns_d(lenbuf),                                 &
                  mat_var(3)%matrix_diagonal_d(n3d) )
         drop_tol=drop_hamiltonian
         Call Pack_Triangle ('hamiltonian',                                            &
                              triangle_hamiltonian_d,                                  &
                              mat_var(3)%packed_columns_d,                             &
                              mat_var(3)%non_zero_columns,                             &
                              mat_var(3)%row_index,                                    &
                              mat_var(3)%number,                                       &
                              mat_var(3)%matrix_diagonal_d)
         DEALLOCATE(triangle_hamiltonian_d,ibuf,h_buf_d,diag_d)
     ELSE IF(ham_type == 'complex') THEN
!
!                  Will not repeat comments as its the same as above but here
!                  for a Hermitian matrix.
!
         ALLOCATE(ibuf(2,max_buf),h_buf_z(max_buf),diag_z(n3d))
         IF (non_orth) THEN
             IF (overlap_to_disk) THEN
                 ALLOCATE(triangle_overlap_z(tri_size))
                 Call Read_and_Write_Packed_Matrices(ibuf,                             &
                                                     h_buf_z,                          &
                                                    'read',                            &
                                                    'overlap',                         &
                                                     non_zero_overlap_elements,        &
                                                     diag_z )
                 Call Fill_Matrix_from_Buffers(triangle_overlap_z,                     &
                                               ibuf,                                   &
                                               h_buf_z,                                &
                                               diag_z,                                 &
                                               non_zero_overlap_elements)
                 write(iout,*) '          *** Overlap Read in.  Non-Zero Elements *** '&
                                              '= ', non_zero_overlap_elements
                 drop_tol=drop_overlap
                 write(iout,*)
                 write(iout,*) '          *** Pack the Overlap ***'
                 ALLOCATE(mat_var(4)%non_zero_columns(n3d),                            &
                          mat_var(4)%row_index(lenbuf),                                &
                          mat_var(4)%packed_columns_z(lenbuf) ,                        &
                          mat_var(4)%matrix_diagonal_z(n3d) )
                 Call Pack_Triangle ('overlap',                                        &
                                     triangle_overlap_z,                               &
                                     mat_var(4)%packed_columns_z,                      &
                                     mat_var(4)%non_zero_columns,                      &
                                     mat_var(4)%row_index,                             &
                                     mat_var(4)%number,                                &
                                     mat_var(4)%matrix_diagonal_z)
                 DEALLOCATE(mat_var(4)%non_zero_columns,                               &
                            mat_var(4)%row_index,                                      &
                            mat_var(4)%packed_columns_z,                               &
                            mat_var(4)%matrix_diagonal_z )
                 DEALLOCATE( triangle_overlap_z )
                 write(iout,*)
             END IF
!
!            The Cholesky factors have already been packed.  Read them in.
!
             write(iout,*) '          *** Reading Packed Upper Cholesky Factor ***'
             Call IOsys('read integer "number_non_zero_upper_elements" from '//        &
                        'packed_matrices',1,mat_var(1)%number,0,' ')
             ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                &
                      mat_var(1)%row_index(mat_var(1)%number),                         &
                      mat_var(1)%packed_columns_z(mat_var(1)%number),                  &
                      mat_var(1)%matrix_diagonal_z(n3d) )
             Call Read_and_Write_Column_Packed_Matrices                                &
                                              (mat_var(1)%packed_columns_z,            &
                                               mat_var(1)%non_zero_columns,            &
                                               mat_var(1)%row_index,                   &
                                               'read',                                 &
                                               'upper',                                &
                                               mat_var(1)%number,                      &
                                               mat_var(1)%matrix_diagonal_z )
             write(iout,*) '          *** Reading Packed Lower Cholesky Factor ***'
             Call IOsys('read integer "number_non_zero_lower_elements" from '//        &
                        'packed_matrices',1,mat_var(2)%number,0,' ')
             ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                &
                      mat_var(2)%row_index(mat_var(2)%number),                         &
                      mat_var(2)%packed_columns_z(mat_var(2)%number) )
             Call Read_and_Write_Column_Packed_Matrices                                &
                                              (mat_var(2)%packed_columns_z,            &
                                               mat_var(2)%non_zero_columns,            &
                                               mat_var(2)%row_index,                   &
                                               'read',                                 &
                                               'lower',                                &
                                               mat_var(2)%number )
         END IF
         ALLOCATE(triangle_hamiltonian_z(tri_size))
         Call Read_and_Write_Packed_Matrices(ibuf,                                     &
                                             h_buf_z,                                  &
                                            'read',                                    &
                                            'hamiltonian',                             &
                                             non_zero_hamiltonian_elements,            & 
                                             diag_z )
         Call Fill_Matrix_from_Buffers(triangle_hamiltonian_z,                         &
                                       ibuf,                                           &
                                       h_buf_z,                                        &
                                       diag_z,                                         &
                                       non_zero_hamiltonian_elements)

         write(iout,*) '          *** Hamiltonian Read in.  Non-Zero Elements *** = ', &
                                                     non_zero_hamiltonian_elements
         write(iout,*)
         write(iout,*) '          *** Pack the Hamiltonian ***'
         ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                    &
                  mat_var(3)%row_index(lenbuf),                                        &
                  mat_var(3)%packed_columns_z(lenbuf),                                 &
                  mat_var(3)%matrix_diagonal_z(n3d) )
         drop_tol=drop_hamiltonian
         Call Pack_Triangle ('hamiltonian',                                            &   
                              triangle_hamiltonian_z,                                  &
                              mat_var(3)%packed_columns_z,                             &
                              mat_var(3)%non_zero_columns,                             &
                              mat_var(3)%row_index,                                    &
                              mat_var(3)%number,                                       &
                              mat_var(3)%matrix_diagonal_z)
         DEALLOCATE(triangle_hamiltonian_z,ibuf,h_buf_z,diag_z)
     END IF
  ELSE IF(output_matrices == 'packed_in_triangular_iosys_format') THEN
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
  ELSE
     Write(iout,1) output_matrices
          Call lnkerr('Bad Output Matrix Keyword')
  END IF
1 FORMAT(/,10X,'Bad Output Matrix Keyword Used = ',a80,/,10x,                         &
               'Allowed Values are: packed_in_standard_iosys_format',/,10x,           &
               '                    packed_in_triangular_iosys_format')
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
  SUBROUTINE Input_Packed_Matrices_d(file_directory,output_matrices)
  IMPLICIT NONE
  REAL*8                                :: gbytes
  REAL*8                                :: dum
  INTEGER                               :: len
  INTEGER                               :: len_1
  INTEGER                               :: lenth
  CHARACTER(LEN=*)                      :: output_matrices
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
                                         .false.)
      stop
  END IF
!
! This is the option to reformat the input matrix in IOsys format and then return 
! to the calling routine.

  write(iout,*) 'Matrices '//output_matrices
  IF(output_matrices == 'packed_in_standard_iosys_format') THEN
!
!                       Standard Two Buffer Random Packing Method
!
     IF(non_orth) THEN      
        READ(50) non_zero_overlap_elements
        write(iout,*) '          *** Pack the Overlap Matrix ***'
        write(iout,*) 'Number Non_Zero Overlap In = ',non_zero_overlap_elements
        ALLOCATE(ibuf(2,non_zero_overlap_elements),                                         &
                 h_buf_d(non_zero_overlap_elements),diag_d(n3d))
        drop_tol=drop_overlap
        Call Fill_Buffers_from_Disk (                                                       &
                                     ibuf,                                                  &
                                     h_buf_d,                                               &
                                     diag_d,                                                &
                                     drop_tol,                                              &
                                     non_zero_overlap_elements,                             &
                                     50)
        write(iout,*) 'Number Non_Zero Overlap Out = ',non_zero_overlap_elements
        write(iout,*) 'Write the Overlap Matrix to Disk'
        IF (overlap_to_disk) THEN
            Call Read_and_Write_Packed_Matrices (                                           &
                                                 ibuf,                                      &
                                                 h_buf_d,                                   &
                                                 'write',                                   &
                                                 'overlap',                                 &
                                                 non_zero_overlap_elements,                 &
                                                 diag_d )
        END IF
!       
!       Cholesky Decomposition and write to DIsk.
!
        write(iout,*)
        write(iout,*)'          *** Form Cholesky Decomposition ***'
        ALLOCATE(upper_d(1:tri_size))
        Call Fill_Matrix_from_Buffers (                                                     &
                                       upper_d,                                             &
                                       ibuf,                                                &
                                       h_buf_d,                                             &
                                       diag_d,                                              &
                                       non_zero_overlap_elements)
        DEALLOCATE(ibuf, h_buf_d )
        call dpptrf('u',n3d,upper_d,info)
        IF(info /= 0) THEN
           write(iout,*) '          *** Singular Overlap Matrix:Quit ***'
           Call lnkerr('Quit.  Singular Overlap Matrix')
        END IF
        IF (full_cholesky_to_disk) THEN
           write(iout,*) '          *** Opening File to Write Full Cholesky Decomposition ***'
           len_1=lenth(cholesky_file_name)
           Call IOsys('open cholesky_decomposition as new',0,0,0,                          &
                       file_directory(1:len)//'/'//cholesky_file_name(1:len_1))
           write(iout,*)'Write Full Cholesky Decomposition to Disk'
           Call IOsys('write real cholesky_factor to cholesky_decomposition',              &
                       tri_size,upper_d,0,' ')
           write(iout,*) 'Closing File Holding Full Cholesky Decomposition'
           Call IOsys('close cholesky_decomposition',0,0,0,' ')
        END IF
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
        ALLOCATE(lower_d(1:tri_size))
        Call Upper_to_Lower(upper_d,lower_d)
        ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                         &
                 mat_var(2)%row_index(lenbuf),                                             &
                 mat_var(2)%packed_columns_d(lenbuf) )
        Call Pack_Triangle ('lower',                                                       &
                             lower_d,                                                      &
                             mat_var(2)%packed_columns_d,                                  &
                             mat_var(2)%non_zero_columns,                                  &
                             mat_var(2)%row_index,                                         &
                             mat_var(2)%number,                                            &
                             mat_var(1)%matrix_diagonal_d)
        DEALLOCATE(mat_var(2)%non_zero_columns,                                            &
                   mat_var(2)%row_index,                                                   &
                   mat_var(2)%packed_columns_d,                                            &
                   mat_var(1)%matrix_diagonal_d )
        DEALLOCATE( upper_d, lower_d )
     END IF
     READ(50) non_zero_hamiltonian_elements
     write(iout,*)
     write(iout,*) '          *** Pack the Hamiltonian Matrix ***'
     write(iout,*) 'Number Non_Zero Hamiltonian In = ',non_zero_hamiltonian_elements
     ALLOCATE(ibuf(2,non_zero_hamiltonian_elements),                                       &
              h_buf_d(non_zero_hamiltonian_elements) )
     drop_tol=drop_hamiltonian
     Call Fill_Buffers_from_Disk (                                                         &
                                  ibuf,                                                    &
                                  h_buf_d,                                                 &
                                  diag_d,                                                  &
                                  drop_tol,                                                &
                                  non_zero_hamiltonian_elements,                           &
                                  50)
     write(iout,*) 'Number Non_Zero Hamiltonian Out = ',non_zero_Hamiltonian_elements
     write(iout,*) '          *** Write Packed Hamiltonian to Disk ***'
     Call Read_and_Write_Packed_Matrices (                                                 &
                                          ibuf,                                            &
                                          h_buf_d,                                         &
                                          'write',                                         &
                                          'hamiltonian',                                   &
                                          non_zero_hamiltonian_elements,                   &
                                          diag_d )
     DEALLOCATE( ibuf, h_buf_d, diag_d )
  ELSE IF(output_matrices == 'packed_in_triangular_iosys_format') THEN
!
!       Triangular Packing by Columns.
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
  ELSE
     Write(iout,1) output_matrices
     Call lnkerr('Bad Output Matrix Keyword')
  END IF
1 FORMAT(/,10X,'Bad Output Matrix Keyword Used = ',a80,/,10x,                              &
               'Allowed Values are: packed_in_standard_iosys_format',/,10x,                &
               '                    packed_in_triangular_iosys_format')
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
  SUBROUTINE Input_Packed_Matrices_z(file_directory,output_matrices)
  IMPLICIT NONE
  REAL*8                                :: gbytes
  COMPLEX*16                            :: dum
  INTEGER                               :: len
  INTEGER                               :: len_1
  INTEGER                               :: lenth
  CHARACTER(LEN=*)                      :: output_matrices
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

  write(iout,*) 'Matrices '//output_matrices
  IF(output_matrices == 'packed_in_standard_iosys_format') THEN
!
!                       Standard Two Buffer Random Packing Method
!
     IF(non_orth) THEN      
        READ(50) non_zero_overlap_elements
        write(iout,*) '          *** Pack the Overlap Matrix ***'
        write(iout,*) 'Number Non_Zero Overlap In = ',non_zero_overlap_elements
        ALLOCATE(ibuf(2,non_zero_overlap_elements),                                         &
                 h_buf_z(non_zero_overlap_elements),diag_z(n3d))
        drop_tol=drop_overlap
        Call Fill_Buffers_from_Disk (                                                       &
                                     ibuf,                                                  &
                                     h_buf_z,                                               &
                                     diag_z,                                                &
                                     drop_tol,                                              &
                                     non_zero_overlap_elements,                             &
                                     50)
        write(iout,*) 'Number Non_Zero Overlap Out = ',non_zero_overlap_elements
        IF (overlap_to_disk) THEN   
           write(iout,*) '          *** Write the Overlap Matrix to Disk ***'
           Call Read_and_Write_Packed_Matrices (                                            &
                                                ibuf,                                       &
                                                h_buf_z,                                    &
                                                'write',                                    &
                                                'overlap',                                  &
                                                non_zero_overlap_elements,                  &
                                                diag_z )
        END IF
!       
!              Cholesky Decomposition and write to DIsk.
!
        write(iout,*)
        write(iout,*)'          *** Form Cholesky Decomposition ***'
        ALLOCATE(upper_z(1:tri_size))
        Call Fill_Matrix_from_Buffers (                                                     &
                                       upper_z,                                             &
                                       ibuf,                                                &
                                       h_buf_z,                                             &
                                       diag_z,                                              &
                                       non_zero_overlap_elements)
        DEALLOCATE(ibuf, h_buf_z)
        call zpptrf('u',n3d,upper_z,info)
        IF(info /= 0) THEN
           write(iout,*) '          *** Singular Overlap Matrix:Quit ***'
           Call lnkerr('Quit.  Singular Overlap Matrix')
        END IF
        IF (full_cholesky_to_disk) THEN
           write(iout,*) '          *** Opening File to Write Full Cholesky Decomposition ***'
           len_1=lenth(cholesky_file_name)
           Call IOsys('open cholesky_decomposition as new',0,0,0,                           &
                       file_directory(1:len)//'/'//cholesky_file_name(1:len_1))
           write(iout,*)'          *** Write Full Cholesky Decomposition to Disk ***'
           Call IOsys('write real cholesky_factor to cholesky_decomposition',               &
                       tri_size,upper_z,0,' ')
           write(iout,*) '          *** Closing File Holding Full Cholesky Decomposition ***'
           Call IOsys('close cholesky_decomposition',0,0,0,' ')
        END IF
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
        ALLOCATE(lower_z(1:tri_size))
        Call Upper_to_Lower(upper_z,lower_z)
        ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                          &
                 mat_var(2)%row_index(lenbuf),                                              &
                 mat_var(2)%packed_columns_z(lenbuf) )
        Call Pack_Triangle ('lower',                                                        &
                             lower_z,                                                       &
                             mat_var(2)%packed_columns_z,                                   &
                             mat_var(2)%non_zero_columns,                                   &
                             mat_var(2)%row_index,                                          &
                             mat_var(2)%number,                                             &
                             mat_var(1)%matrix_diagonal_z)
        DEALLOCATE(mat_var(2)%non_zero_columns,                                             &
                   mat_var(2)%row_index,                                                    &
                   mat_var(2)%packed_columns_z,                                             &
                   mat_var(1)%matrix_diagonal_z )
        DEALLOCATE( upper_z, lower_z)
     END IF
     READ(50) non_zero_hamiltonian_elements
     write(iout,*)
     write(iout,*) '          *** Pack the Hamiltonian Matrix ***'
     write(iout,*) 'Number Non_Zero Hamiltonian In = ',non_zero_hamiltonian_elements
     ALLOCATE(ibuf(2,non_zero_hamiltonian_elements),                                        &
              h_buf_z(non_zero_hamiltonian_elements) )
     drop_tol=drop_hamiltonian
     Call Fill_Buffers_from_Disk (                                                          &
                                  ibuf,                                                     &
                                  h_buf_z,                                                  &
                                  diag_z,                                                   &
                                  drop_tol,                                                 &
                                  non_zero_hamiltonian_elements,                            &
                                  50)
     write(iout,*) 'Number Non_Zero Hamiltonian Out = ',non_zero_hamiltonian_elements
     write(iout,*) '          *** Write Packed Hamiltonian to Disk ***'
     Call Read_and_Write_Packed_Matrices (                                                  &
                                          ibuf,                                             &
                                          h_buf_z,                                          &
                                          'write',                                          &
                                          'hamiltonian',                                    &
                                          non_zero_hamiltonian_elements,                    &
                                          diag_z )
     DEALLOCATE( ibuf, h_buf_z, diag_z )
  ELSE IF(output_matrices == 'packed_in_triangular_iosys_format') THEN
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
  ELSE
     Write(iout,1) output_matrices
     Call lnkerr('Bad Output Matrix Keyword')
  END IF
1 FORMAT(/,10X,'Bad Output Matrix Keyword Used = ',a80,/,10x,                              &
               'Allowed Values are: packed_in_standard_iosys_format',/,10x,                &
               '                    packed_in_triangular_iosys_format')
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
