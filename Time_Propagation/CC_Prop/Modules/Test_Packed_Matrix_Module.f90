!***********************************************************************
                           MODULE Packed_Matrix_Module
!
                           USE Global_Time_Propagation_Module
                           USE Pack_Hamiltonian_Module
                           USE Preconditioner_Module
!
                           IMPLICIT NONE
            CHARACTER(LEN=1600)                  :: data_card
            CHARACTER(LEN=80)                    :: pass_data
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
      len=lenth(file_directory)
      len_1=lenth(ham_file)
      write(iout,*) '*** Opening and Processing Input Matrix File ***'
      OPEN(UNIT=50,FILE=file_directory(1:len)//'/'//ham_file(1:len_1),               &
           ACCESS='sequential', FORM='unformatted',IOSTAT=IOSTAT,STATUS='old')
      READ(50)
      write(iout,*) '    *** Opening File to Write Packed Matrices ***'
      Call IOsys('open packed_matrices as new',0,0,0,                                &
                  file_directory(1:len)//'/'//'Packed_Matrices')
      IF (ham_type == 'real') THEN
          Call Input_Packed_Matrices_d(file_directory,output_matrices)
      ELSE IF(ham_type == 'complex') THEN
          Call Input_Packed_Matrices_z(file_directory,output_matrices)
      END IF
      write(iout,*)
      write(iout,*) '*** Closing Input Matrix File ***'
      CLOSE(50)
      write(iout,*) '*** Ending Writing of Packed Matrices ***'
      Call IOsys ('close packed_matrices',0,0,0,' ')
      IF( reformat_input_matrix_only) THEN
!
!         We can opt to quit now if we wish or to re-read in the matrices for
!         a full calculation.
!
          write(iout,*) '*** Quitting After Reformat ***'
          stop
      ELSE
          write(iout,*)
          write(iout,*) '*** Reopening File to Read Packed Matrices ***'
          Call IOsys('open packed_matrices as old',0,0,0,                            &
                      file_directory(1:len)//'/'//'Packed_Matrices')
          Call Input_Packed_IOsys_Matrix(ham_type,file_directory,output_matrices)
          write(iout,*) '***Closing File Packed Matrices ***'
          Call IOsys('close packed_matrices',0,0,0,' ')
      END IF
  ELSE IF(input_matrices =='packed_in_iosys_format') THEN
!
!         Take the already IOsys formatted files and read in the matrices.
!
          len=lenth(file_directory)
          write(iout,*) '*** Opening File to Read Packed Matrices ***'
          Call IOsys('open packed_matrices as old',0,0,0,                            &
                      file_directory(1:len)//'/'//'Packed_Matrices')
          Call Input_Packed_IOsys_Matrix(ham_type,file_directory,output_matrices)
          write(iout,*) '*** Closing File Packed Matrices ***'
          Call IOsys('close packed_matrices',0,0,0,' ')
  END IF
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
     write(iout,*) '***Reading in Packed Matrices in Standard Two Buffer Form ***'
     Call IOsys('read integer number_non_zero_hamiltonian_elements from '//            &
                'packed_matrices',1,non_zero_hamiltonian_elements,0,' ') 
     IF (non_orth) THEN
         Call IOsys('read integer number_non_zero_overlap_elements from '//            &
                    'packed_matrices',1,non_zero_overlap_elements,0,' ')
         Call IOsys('read integer number_non_zero_cholesky_elements from '//           &
                    'packed_matrices',1,non_zero_cholesky_elements,0,' ') 
     END IF
     max_buf = max(non_zero_overlap_elements,non_zero_hamiltonian_elements,            &
                   non_zero_cholesky_elements)
!
!          In this section we read in the already packed matrices and then
!          fill up the full matrices with the elements.
!          
     IF (ham_type == 'real') THEN
         ALLOCATE(ibuf(2,max_buf),h_buf_d(max_buf),diag_d(n3d))
         IF (non_orth) THEN
             ALLOCATE(triangle_overlap_d(tri_size))
!
!            Read in the packed matrices into the buffers.
!
             Call Read_and_Write_Packed_Matrices(ibuf,                                 &
                                                 h_buf_d,                              &
                                                 diag_d,                               &
                                                 non_zero_overlap_elements,            &
                                                'overlap',                             &
                                                'read')
!
!            Fill the full matrix from the buffers.
!
             triangle_overlap_d(1:tri_size) = 0.d0 
             Call Fill_Matrix_from_Buffers(triangle_overlap_d,                         &
                                           ibuf,                                       &
                                           h_buf_d,                                    &
                                           diag_d,                                     &
                                           non_zero_overlap_elements)
             write(iout,*) '   *** Overlap Read in.  Non-Zero Elements *** = ',        &
                                                     non_zero_overlap_elements
!
!            Form the column packed representation and write it out to disk.
!
             drop_tol=drop_overlap
             write(iout,*)
             write(iout,*) '   *** Pack the Overlap ***'
             ALLOCATE(mat_var(4)%non_zero_columns(n3d),                                &
                      mat_var(4)%row_index(lenbuf),                                    &
                      mat_var(4)%packed_columns_d(lenbuf) ,                            &
                      mat_var(4)%matrix_diagonal_d(n3d) )
             Call Pack_Triangle_d(triangle_overlap_d,'overlap')
             DEALLOCATE(mat_var(4)%non_zero_columns,                                   &
                        mat_var(4)%row_index,                                          &
                        mat_var(4)%packed_columns_d,                                   &
                        mat_var(4)%matrix_diagonal_d )
!
!            Do similar business for Cholesky factors
!
             write(iout,*)
             write(iout,*) '   *** Pack the Cholesky Factors ***'
             ALLOCATE(upper_d(tri_size))
             Call Read_and_Write_Packed_Matrices(ibuf,                                 &
                                                 h_buf_d,                              &
                                                 diag_d,                               &
                                                 non_zero_cholesky_elements,           &
                                                'cholesky',                            &
                                                'read')
             upper_d(1:tri_size) = 0.d0 
             Call Fill_Matrix_from_Buffers(upper_d,ibuf,h_buf_d,diag_d,                &
                              non_zero_cholesky_elements)
             write(iout,*) '   *** Cholesky Read in.  Non-Zero Elements *** = ',       &
                                                      non_zero_cholesky_elements
             drop_tol=drop_cholesky
             ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                &
                      mat_var(1)%row_index(lenbuf),                                    &
                      mat_var(1)%packed_columns_d(lenbuf) ,                            &
                      mat_var(1)%matrix_diagonal_d(n3d) )
             Call Pack_Triangle(upper_d,'upper')
             Call Upper_to_Lower(upper_d,triangle_overlap_d)
             DEALLOCATE(upper_d)
             ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                &
                      mat_var(2)%row_index(lenbuf),                                    &
                      mat_var(2)%packed_columns_d(lenbuf) )
             Call Pack_Triangle(triangle_overlap_d,'lower')
             DEALLOCATE(triangle_overlap_d)
         END IF
!
!        Do similarly for Hamiltonian
!
         ALLOCATE(triangle_hamiltonian_d(tri_size))
         triangle_hamiltonian_d(1:tri_size) = 0.d0 
         Call Read_and_Write_Packed_Matrices(ibuf,                                     &
                                             h_buf_d,                                  &
                                             diag_d,                                   &
                                             non_zero_hamiltonian_elements,            & 
                                            'hamiltonian',                             &
                                            'read')
         Call Fill_Matrix_from_Buffers(triangle_hamiltonian_d,                         &
                                       ibuf,                                           &
                                       h_buf_d,                                        &
                                       diag_d,                                         &
                                       non_zero_hamiltonian_elements)
         write(iout,*) '   *** Hamiltonian Read in.  Non-Zero Elements *** = ',        &
                                                     non_zero_hamiltonian_elements
         write(iout,*)
         write(iout,*) '   *** Pack the Hamiltonian ***'
         ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                    &
                  mat_var(3)%row_index(lenbuf),                                        &
                  mat_var(3)%packed_columns_d(lenbuf),                                 &
                  mat_var(3)%matrix_diagonal_d(n3d) )
         drop_tol=drop_hamiltonian
         Call Pack_Triangle(triangle_hamiltonian_d,'hamiltonian')
         DEALLOCATE(triangle_hamiltonian_d,ibuf,h_buf_d,diag_d)
     ELSE IF(ham_type == 'complex') THEN
!
!                  Will not repeat comments as its the same as above but here
!                  for a Hermitian matrix.
!
         ALLOCATE(ibuf(2,max_buf),h_buf_z(max_buf),diag_z(n3d))
         IF (non_orth) THEN
             ALLOCATE(triangle_overlap_z(tri_size))
             Call Read_and_Write_Packed_Matrices(ibuf,                                 &
                                                 h_buf_z,                              &
                                                 diag_z,                               &
                                                 non_zero_overlap_elements,            &
                                                'overlap',                             &
                                                'read')
             triangle_overlap_z(1:tri_size) = 0.d0 
             Call Fill_Matrix_from_Buffers(triangle_overlap_z,                         &
                                           ibuf,                                       &
                                           h_buf_z,                                    &
                                           diag_z,                                     &
                                           non_zero_overlap_elements)
             write(iout,*) '   *** Overlap Read in.  Non-Zero Elements *** = ',        &
                                                     non_zero_overlap_elements
             drop_tol=drop_overlap
             write(iout,*)
             write(iout,*) '   *** Pack the Overlap ***'
             ALLOCATE(mat_var(4)%non_zero_columns(n3d),                                &
                      mat_var(4)%row_index(lenbuf),                                    &
                      mat_var(4)%packed_columns_z(lenbuf) ,                            &
                      mat_var(4)%matrix_diagonal_z(n3d) )
             Call Pack_Triangle_z(triangle_overlap_z,'overlap')
             DEALLOCATE(mat_var(4)%non_zero_columns,                                   &
                        mat_var(4)%row_index,                                          &
                        mat_var(4)%packed_columns_z,                                   &
                        mat_var(4)%matrix_diagonal_z )
             write(iout,*)
             write(iout,*) '   *** Pack the Cholesky Factors ***'
             ALLOCATE(upper_z(tri_size))
             Call Read_and_Write_Packed_Matrices(ibuf,                                 &
                                                 h_buf_z,                              &
                                                 diag_z,                               &
                                                 non_zero_cholesky_elements,           &
                                                'cholesky',                            &
                                                'read')
             upper_z(1:tri_size) = 0.d0 
             Call Fill_Matrix_from_Buffers(upper_z,                                    &
                                           ibuf,                                       &
                                           h_buf_z,                                    &
                                           diag_z,                                     &
                                           non_zero_cholesky_elements)
             write(iout,*) '   *** Cholesky Read in.  Non-Zero Elements *** = ',       &
                                                      non_zero_cholesky_elements
             drop_tol=drop_cholesky
             ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                &
                      mat_var(1)%row_index(lenbuf),                                    &
                      mat_var(1)%packed_columns_z(lenbuf) ,                            &
                      mat_var(1)%matrix_diagonal_z(n3d) )
             Call Pack_Triangle(upper_z,'upper')
             Call Upper_to_Lower(upper_z,triangle_overlap_z)
             DEALLOCATE(upper_z)
             ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                &
                      mat_var(2)%row_index(lenbuf),                                    &
                      mat_var(2)%packed_columns_z(lenbuf) )
             Call Pack_Triangle(triangle_overlap_z,'lower')
             DEALLOCATE(triangle_overlap_z)
         END IF
         ALLOCATE(triangle_hamiltonian_z(tri_size))
         triangle_hamiltonian_z(1:tri_size) = 0.d0 
         Call Read_and_Write_Packed_Matrices(ibuf,                                     &
                                             h_buf_z,                                  &
                                             diag_z,                                   &
                                             non_zero_hamiltonian_elements,            & 
                                            'hamiltonian',                             &
                                            'read')
         Call Fill_Matrix_from_Buffers(triangle_hamiltonian_z,                         &
                                       ibuf,                                           &
                                       h_buf_z,                                        &
                                       diag_z,                                         &
                                       non_zero_hamiltonian_elements)

         write(iout,*) '   *** Hamiltonian Read in.  Non-Zero Elements *** = ',        &
                                                     non_zero_hamiltonian_elements
         write(iout,*)
         write(iout,*) '   *** Pack the Hamiltonian ***'
         ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                    &
                  mat_var(3)%row_index(lenbuf),                                        &
                  mat_var(3)%packed_columns_z(lenbuf),                                 &
                  mat_var(3)%matrix_diagonal_z(n3d) )
         drop_tol=drop_hamiltonian
         Call Pack_Triangle(triangle_hamiltonian_z,'hamiltonian')
         DEALLOCATE(triangle_hamiltonian_z,ibuf,h_buf_z,diag_z)
     END IF
  ELSE
     write(iout,*) '*** Reading in Packed Matrices in Triangle Form ***'
!
!    The matrices here are already packed in column form and only need to be read in
!
     IF (ham_type == 'real') THEN
         IF (non_orth) THEN
             write(iout,*) '   *** Reading Packed Upper Cholesky Factor ***'
             ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                &
                      mat_var(1)%row_index(lenbuf),                                    &
                      mat_var(1)%packed_columns_d(lenbuf) ,                            &
                      mat_var(1)%matrix_diagonal_d(n3d) )
             Call Read_and_Write_Column_Packed_Matrices                                &
                                              (mat_var(1)%packed_columns_d,            &
                                               mat_var(1)%non_zero_columns,            &
                                               mat_var(1)%row_index,                   &
                                               mat_var(1)%number,                      &
                                               'upper',                                &
                                               'read',                                 &
                                               mat_var(1)%matrix_diagonal_d )
             write(iout,*) '   *** Reading Packed Lower Cholesky Factor ***'
             ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                &
                      mat_var(2)%row_index(lenbuf),                                    &
                      mat_var(2)%packed_columns_d(lenbuf) )
             Call Read_and_Write_Column_Packed_Matrices                                &
                                              (mat_var(2)%packed_columns_d,            &
                                               mat_var(2)%non_zero_columns,            &
                                               mat_var(2)%row_index,                   &
                                               mat_var(2)%number,                      &
                                               'lower',                                &
                                               'read')
         END IF
         write(iout,*) '   *** Reading Packed Hamiltonian ***'
         ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                    &
                  mat_var(3)%row_index(lenbuf),                                        &
                  mat_var(3)%packed_columns_d(lenbuf),                                 &
                  mat_var(3)%matrix_diagonal_d(n3d) )
         Call Read_and_Write_Column_Packed_Matrices                                    &
                                          (mat_var(3)%packed_columns_d,                &
                                           mat_var(3)%non_zero_columns,                &
                                           mat_var(3)%row_index,                       &
                                           mat_var(3)%number,                          &
                                          'hamiltonian',                               &
                                          'read',                                      &
                                           mat_var(3)%matrix_diagonal_d )
     ELSE IF(ham_type == 'complex') THEN
        IF (non_orth) THEN
            write(iout,*) 'Reading Packed Upper Cholesky Factor'
            ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                 &
                     mat_var(1)%row_index(lenbuf),                                     &
                     mat_var(1)%packed_columns_z(lenbuf) ,                             &
                     mat_var(1)%matrix_diagonal_z(n3d) )
            Call Read_and_Write_Column_Packed_Matrices                                 &
                                             (mat_var(1)%packed_columns_z,             &
                                              mat_var(1)%non_zero_columns,             &
                                              mat_var(1)%row_index,                    &
                                              mat_var(1)%number,                       &
                                              'upper',                                 &
                                              'read',                                  &
                                              mat_var(1)%matrix_diagonal_z )
            write(iout,*) '   *** Reading Packed Lower Cholesky Factor ***'
            ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                 &
                     mat_var(2)%row_index(lenbuf),                                     &
                     mat_var(2)%packed_columns_z(lenbuf) )
            Call Read_and_Write_Column_Packed_Matrices                                 &
                                             (mat_var(2)%packed_columns_z,             &
                                              mat_var(2)%non_zero_columns,             &
                                              mat_var(2)%row_index,                    &
                                              mat_var(2)%number,                       &
                                              'lower',                                 &
                                              'read')
        END IF
        write(iout,*) '   *** Reading Packed Hamiltonian ***'
        ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                     &
                 mat_var(3)%row_index(lenbuf),                                         &
                 mat_var(3)%packed_columns_z(lenbuf),                                  &
                 mat_var(3)%matrix_diagonal_z(n3d) )
        Call Read_and_Write_Column_Packed_Matrices                                     &
                                         (mat_var(3)%packed_columns_z,                 &
                                          mat_var(3)%non_zero_columns,                 &
                                          mat_var(3)%row_index,                        &
                                          mat_var(3)%number,                           &
                                         'hamiltonian',                                &
                                         'read',                                       &
                                         mat_var(3)%matrix_diagonal_z )
     END IF
  END IF
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
  REAL*8                                :: matrix_el
!  INTEGER                               :: i
!  INTEGER                               :: j
  REAL*8                                :: gbytes
  REAL*8                                :: dum
  INTEGER                               :: len
  INTEGER                               :: lenth
  CHARACTER(LEN=*)                      :: output_matrices
  CHARACTER(LEN=*)                      :: file_directory
  gbytes = tri_size*8*3
!
! This option is just used for test purposes and if this option is set, we exit 
! immediately after the diagonalization.
!
  len = lenth(file_directory)
  IF (diagonalize_only) THEN
      write(iout,*) '*** Diagonalizing Big Matrix ***'
      IF(non_orth) THEN
         lwork = 10 * n3d
         ALLOCATE( triangle_hamiltonian_d(1:tri_size), triangle_overlap_d(1:tri_size),        &
                                          eig(1:n3d),  work_d(1:lwork) )
         Call Fill_Matrix_from_Disk(triangle_overlap_d,                                       &
                                    non_zero_overlap_elements,                                &
                                    50)
         Call Fill_Matrix_from_Disk(triangle_hamiltonian_d,                                   &
                                    non_zero_hamiltonian_elements,                            &
                                    50)
         call dspgv(1,'n','u',n3d,triangle_hamiltonian_d,                                     &
                                  triangle_overlap_d,                                         &
                                  eig,dum,n3d,work_d,lwork,info)
      ELSE
         lwork=3*n3d 
         ALLOCATE( triangle_hamiltonian_d(1:tri_size), eig(1:n3d), work_d(1:lwork) )
         Call Fill_Matrix_from_Disk(triangle_hamiltonian_d,                                   &
                                    non_zero_hamiltonian_elements,                            &
                                    50)
         call dspev('n','u',n3d,triangle_hamiltonian_d,eig,dum,n3d,work_d,info)
      END IF
      write(iout,*) '*** Done Diagonalizing Big Matrix: Quit ***'
      local_title='eigenvalues'
      call prntfmn(local_title,eig,n3d,1,n3d,1,iout,'e')
      Call IOsys('write real eigenvalues to packed_matrices',n3d,eig,0,' ')
      Call IOsys('close packed_matrices',0,0,0,' ')
      stop
  END IF
!
! This is the option to reformat the input matrix in IOsys format and then 
! return to the calling routine.

  write(iout,*) 'Matrices '//output_matrices
  IF(output_matrices == 'packed_in_standard_iosys_format') THEN
!
!            Here we use the Standard Two Buffer Random Packing Method
!
     IF(non_orth) THEN      
        READ(50) non_zero_overlap_elements
        write(iout,*) '   *** Pack the Overlap Matrix'
        write(iout,*) '       *** Number Non_Zero Overlap In *** = ',                         &
                                  non_zero_overlap_elements
!
!                 Fill the buffers from the input diak
!
        ALLOCATE(ibuf(2,non_zero_overlap_elements),                                           &
                 h_buf_d(non_zero_overlap_elements),diag_d(n3d))
!
!                 Overlap first
!
        drop_tol=drop_overlap
        Call Fill_Buffers_from_Disk(ibuf,                                                     &
                                    h_buf_d,                                                  &
                                    diag_d,                                                   &
                                    drop_tol,                                                 &
                                    non_zero_overlap_elements,                                &
                                    50)
        write(iout,*) '       *** Number Non_Zero Overlap Out *** = ',                        &
                                 non_zero_overlap_elements
        write(iout,*) '   *** Write the Overlap Matrix to Disk ***'
!
!                  Write the buffers to disk.
! 
        Call Read_and_Write_Packed_Matrices(ibuf,                                             &
                                            h_buf_d,                                          &
                                            diag_d,                                           &
                                            non_zero_overlap_elements,                        &
                                            'overlap',                                        &
                                            'write')
!       
!              Perform the Cholesky Decomposition and write to Disk in both full
!              and column packed form..
!
        write(iout,*)
        write(iout,*)'   *** Form Cholesky Decomposition ***'
        ALLOCATE(upper_d(1:tri_size))
        Call Fill_Matrix_from_Buffers(upper_d,                                                &
                                      ibuf,                                                   &
                                      h_buf_d,                                                &
                                      diag_d,                                                 &
                                      non_zero_overlap_elements)
        call dpptrf('u',n3d,upper_d,info)
        IF(info /= 0) THEN
           write(iout,*) '   *** Singular Overlap Matrix:Quit ***'
           Call lnkerr('Quit.  Singular Overlap Matrix')
        END IF
        write(iout,*) '   *** Opening File to Write Full Cholesky Decomposition ***'
        Call IOsys('open cholesky_decomposition as new',0,0,0,                                &
                    file_directory(1:len)//'/'//'Cholesky_Decomposition')
        write(iout,*)'      *** Write Full Cholesky Decomposition to Disk ***'
        Call IOsys('write real cholesky_factor to cholesky_decomposition',                    &
                    tri_size,upper_d,0,' ')
        write(iout,*) '   *** Closing File Holding Full Cholesky Decomposition ***'
        Call IOsys('close cholesky_decomposition',0,0,0,' ')
        drop_tol=drop_cholesky
        ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                            &
                 mat_var(1)%row_index(lenbuf),                                                &
                 mat_var(1)%packed_columns_d(lenbuf) ,                                        &
                 mat_var(1)%matrix_diagonal_d(n3d) )
        Call Pack_Triangle(upper_d,'upper')
        DEALLOCATE(mat_var(1)%non_zero_columns,                                               &
                   mat_var(1)%row_index,                                                      &
                   mat_var(1)%packed_columns_d )
        ALLOCATE(lower_d(1:tri_size))
        Call Upper_to_Lower(upper_d,lower_d)                                      
        ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                            &
                 mat_var(2)%row_index(lenbuf),                                                &
                 mat_var(2)%packed_columns_d(lenbuf) )
        Call Pack_Triangle(lower_d,'lower')
        DEALLOCATE(mat_var(2)%non_zero_columns,                                               &
                   mat_var(2)%row_index,                                                      &
                   mat_var(2)%packed_columns_d,                                               &
                   mat_var(1)%matrix_diagonal_d )
        DEALLOCATE(ibuf, h_buf_d, lower_d)
        non_zero_cholesky_elements = tri_size
        ALLOCATE(ibuf(2,non_zero_cholesky_elements),                                          &
                 h_buf_d(non_zero_cholesky_elements))
        Call Fill_Buffers_from_Matrix(upper_d,                                                &
                                      ibuf,                                                   &
                                      h_buf_d,                                                &
                                      diag_d,                                                 &
                                      drop_tol,                                               &
                                      non_zero_cholesky_elements)
        DEALLOCATE(upper_d)
        write(iout,*)
        write(iout,*) '   *** Pack the Cholesky Factors ***'
        write(iout,*) '      *** Number Non_Zero Cholesky Out *** = ',                        &
                                 non_zero_cholesky_elements 
        write(iout,*)'    *** Write Packed Cholesky Factors to Disk ***'
        Call Read_and_Write_Packed_Matrices(ibuf,                                             &
                                            h_buf_d,                                          &
                                            diag_d,                                           &
                                            non_zero_cholesky_elements,                       &
                                           'cholesky',                                        &
                                           'write')
        DEALLOCATE(ibuf, h_buf_d, diag_d )
     END IF
!
!               Now finish with the Hamiltonian matrix elements.
!
     READ(50) non_zero_hamiltonian_elements
     write(iout,*)
     write(iout,*) '   *** Pack the Hamiltonian Matrix ***'
     write(iout,*) '      *** Number Non_Zero Hamiltonian In *** = ',                         &
                              non_zero_hamiltonian_elements
     drop_tol=drop_hamiltonian
     ALLOCATE(ibuf(2,non_zero_hamiltonian_elements),                                          &
              h_buf_d(non_zero_hamiltonian_elements),diag_d(n3d))
     Call Fill_Buffers_from_Disk(ibuf,                                                        &
                                 h_buf_d,                                                     &
                                 diag_d,                                                      &
                                 drop_tol,                                                    &
                                 non_zero_hamiltonian_elements,                               &
                                 50)
     write(iout,*) '      *** Number Non_Zero Hamiltonian Out *** = ',                        &
                              non_zero_hamiltonian_elements
     write(iout,*) '   *** Write Packed Hamiltonian to Disk ***'
     Call Read_and_Write_Packed_Matrices(ibuf,                                                &
                                         h_buf_d,                                             &
                                         diag_d,                                              &
                                         non_zero_hamiltonian_elements,                       &
                                         'hamiltonian',                                       &
                                         'write')
     DEALLOCATE( ibuf, h_buf_d, diag_d )
  ELSE
!
!              This section does the Column packing directly.
!
     IF(non_orth) THEN      
        ALLOCATE(triangle_overlap_d(tri_size))
        Call Fill_Matrix_from_Disk(triangle_overlap_d,                                        &
                                   non_zero_overlap_elements,                                 &
                                   50)
        write(iout,*) '      *** Number Non_Zero Overlap In *** = ',                          &
                                 non_zero_overlap_elements
        drop_tol=drop_overlap
        write(iout,*)
        write(iout,*) '   *** Pack the Overlap ***'
        ALLOCATE(mat_var(4)%non_zero_columns(n3d),                                           &
                 mat_var(4)%row_index(lenbuf),                                               &
                 mat_var(4)%packed_columns_d(lenbuf) ,                                       &
                 mat_var(4)%matrix_diagonal_d(n3d) )
        print_packed_matrices=.true.
        print_internal_matrices=.true.
!        write(iout,*) 'Write Overlap'
!        write(iout,*) triangle_overlap_d
        Call Pack_Triangle_d(triangle_overlap_d,'overlap')
!       write(iout,*) 'Write Overlap'
!        write(iout,*) triangle_overlap_d
        Call Read_and_Write_Column_Packed_Matrices                                           &
                                         (mat_var(4)%packed_columns_d,                       &
                                          mat_var(4)%non_zero_columns,                       &
                                          mat_var(4)%row_index,                              &
                                          mat_var(4)%number,                                 &
                                          'overlap',                                         &
                                          'read',                                            &
                                          mat_var(4)%matrix_diagonal_d )
        Call Fill_Upper_Triangular_Matrix_from_Column_Buffers_d (                            &
                                           triangle_overlap_d,                               &
                                           mat_var(4)%matrix_diagonal_d,                     &
                                           mat_var(4)%packed_columns_d,                      &
                                           mat_var(4)%non_zero_columns,                      &
                                           mat_var(4)%row_index)
     ALLOCATE(triangle_hamiltonian_d(tri_size))
     Call Fill_Matrix_from_Disk(triangle_hamiltonian_d,                                      &
                                non_zero_hamiltonian_elements,                               &
                                50)
     write(iout,*) '      *** Number Non_Zero Hamiltonian In *** = ',                        &
                              non_zero_hamiltonian_elements
     write(iout,*) 'Write Hamiltonian'
     write(iout,*) triangle_hamiltonian_d
     write(iout,*)
     write(iout,*) '   *** Pack the Hamiltonian ***'
     ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                             &
              mat_var(3)%row_index(lenbuf),                                                 &
              mat_var(3)%packed_columns_d(lenbuf),                                          &
              mat_var(3)%matrix_diagonal_d(n3d) )
     drop_tol=drop_hamiltonian
     Call Pack_Triangle(triangle_hamiltonian_d,'hamiltonian')
     Call Read_and_Write_Column_Packed_Matrices                                             &
                                         (mat_var(3)%packed_columns_d,                      &
                                          mat_var(3)%non_zero_columns,                      &
                                          mat_var(3)%row_index,                             &
                                          mat_var(3)%number,                                &
                                          'hamiltonian',                                    &
                                          'read',                                           &
                                          mat_var(3)%matrix_diagonal_d )
     Call Fill_Upper_Triangular_Matrix_from_Column_Buffers_d (                              &
                                           triangle_hamiltonian_d,                          &
                                           mat_var(3)%matrix_diagonal_d,                    &
                                           mat_var(3)%packed_columns_d,                     &
                                           mat_var(3)%non_zero_columns,                     &
                                           mat_var(3)%row_index)
         lwork=10*n3d
         ALLOCATE( eig(1:n3d), work_d(1:lwork) )
         call dspgv(1,'n','u',n3d,triangle_hamiltonian_d,                                   &
                                  triangle_overlap_d,                                       &
                                  eig,dum,n3d,work_d,lwork,info)
         local_title='eigenvalues'
         call prntfmn(local_title,eig,n3d,1,n3d,1,iout,'e')
         stop
     DEALLOCATE(mat_var(3)%non_zero_columns,                                                &
                mat_var(3)%row_index,                                                       &
                mat_var(3)%packed_columns_d,                                                &
                mat_var(3)%matrix_diagonal_d )
     DEALLOCATE(triangle_hamiltonian_d)
         stop




        DEALLOCATE(mat_var(4)%non_zero_columns,                                             &
                   mat_var(4)%row_index,                                                    &
                   mat_var(4)%packed_columns_d,                                             &
                   mat_var(4)%matrix_diagonal_d )
!
!                 Perform Cholesky decomposition
!
        ALLOCATE(upper_d(tri_size))
        upper_d(1:tri_size) = triangle_overlap_d(1:tri_size)
        write(iout,*)
        write(iout,*)'   *** Form Cholesky Decomposition ***'
        call dpptrf('u',n3d,upper_d,info)
        IF(info /= 0) THEN
           write(iout,*) '   *** Singular Overlap Matrix:Quit ***'
           Call lnkerr('Quit.  Singular Overlap Matrix')
        END IF
        write(iout,*)
        write(iout,*) '   *** Pack the Cholesky Factors ***'
        drop_tol=drop_cholesky
        ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                          &
                 mat_var(1)%row_index(lenbuf),                                              &
                 mat_var(1)%packed_columns_d(lenbuf) ,                                      &
                 mat_var(1)%matrix_diagonal_d(n3d) )
        Call Pack_Triangle(upper_d,'upper')
        DEALLOCATE(mat_var(1)%non_zero_columns,                                             &
                   mat_var(1)%row_index,                                                    &
                   mat_var(1)%packed_columns_d )
        Call Upper_to_Lower(upper_d,triangle_overlap_d)
        DEALLOCATE(upper_d)
        ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                          &
                 mat_var(2)%row_index(lenbuf),                                              &
                 mat_var(2)%packed_columns_d(lenbuf) )
        Call Pack_Triangle(triangle_overlap_d,'lower')
        DEALLOCATE(mat_var(2)%non_zero_columns,                                             &
                   mat_var(2)%row_index,                                                    &
                   mat_var(2)%packed_columns_d,                                             &
                   mat_var(1)%matrix_diagonal_d )
        DEALLOCATE(triangle_overlap_d)
!
!                    Test Section
        ALLOCATE(mat_var(4)%non_zero_columns(n3d),                                          &
                 mat_var(4)%row_index(lenbuf),                                              &
                 mat_var(4)%packed_columns_d(lenbuf) ,                                      &
                 mat_var(4)%matrix_diagonal_d(n3d) )
        ALLOCATE(triangle_overlap_d(1:tri_size))
        Call Read_and_Write_Column_Packed_Matrices                                          &
                                         (mat_var(4)%packed_columns_d,                      &
                                          mat_var(4)%non_zero_columns,                      &
                                          mat_var(4)%row_index,                             &
                                          mat_var(4)%number,                                &
                                          'overlap',                                        &
                                          'read',                                           &
                                          mat_var(4)%matrix_diagonal_d )
        Call Fill_Upper_Triangular_Matrix_from_Column_Buffers_d (                           &
                                           triangle_overlap_d,                              &
                                           mat_var(4)%matrix_diagonal_d,                    &
                                           mat_var(4)%packed_columns_d,                     &
                                           mat_var(4)%non_zero_columns,                     &
                                           mat_var(4)%row_index)
         stop
     END IF
!
!                Finish with the Hamiltonian matrix.
!
     ALLOCATE(triangle_hamiltonian_d(tri_size))
     Call Fill_Matrix_from_Disk(triangle_hamiltonian_d,                                     &
                                non_zero_hamiltonian_elements,                              &
                                50)
     write(iout,*) '      *** Number Non_Zero Hamiltonian In *** = ',                       &
                              non_zero_hamiltonian_elements
     write(iout,*)
     write(iout,*) '   *** Pack the Hamiltonian ***'
     ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                             &
              mat_var(3)%row_index(lenbuf),                                                 &
              mat_var(3)%packed_columns_d(lenbuf),                                          &
              mat_var(3)%matrix_diagonal_d(n3d) )
     drop_tol=drop_hamiltonian
     Call Pack_Triangle(triangle_hamiltonian_d,'hamiltonian')
     DEALLOCATE(mat_var(3)%non_zero_columns,                                                &
                mat_var(3)%row_index,                                                       &
                mat_var(3)%packed_columns_d,                                                &
                mat_var(3)%matrix_diagonal_d )
     DEALLOCATE(triangle_hamiltonian_d)
     stop
!
!            Test Section
!
        ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                          &
                 mat_var(3)%row_index(lenbuf),                                              &
                 mat_var(3)%packed_columns_d(lenbuf) ,                                      &
                 mat_var(3)%matrix_diagonal_d(n3d) )
        ALLOCATE(triangle_hamiltonian_d(1:tri_size))
        Call Read_and_Write_Column_Packed_Matrices                                          &
                                         (mat_var(3)%packed_columns_d,                      &
                                          mat_var(3)%non_zero_columns,                      &
                                          mat_var(3)%row_index,                             &
                                          mat_var(3)%number,                                &
                                          'hamiltonian',                                    &
                                          'read',                                           &
                                          mat_var(3)%matrix_diagonal_d )
        Call Fill_Upper_Triangular_Matrix_from_Column_Buffers_d (                           &
                                           triangle_hamiltonian_d,                          &
                                           mat_var(3)%matrix_diagonal_d,                    &
                                           mat_var(3)%packed_columns_d,                     &
                                           mat_var(3)%non_zero_columns,                     &
                                           mat_var(3)%row_index)
         write(iout,*) triangle_hamiltonian_d
         lwork=10*n3d
         ALLOCATE( eig(1:n3d), work_d(1:lwork) )
         call dspgv(1,'n','u',n3d,triangle_hamiltonian_d,                                   &
                                  triangle_overlap_d,                                       &
                                  eig,dum,n3d,work_d,lwork,info)
         local_title='eigenvalues'
         call prntfmn(local_title,eig,n3d,1,n3d,1,iout,'e')
         DEALLOCATE(triangle_overlap_d)
         stop
  END IF
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
!***purpose            Read in the lower triangle of real symmetric matrix
!***
!***references
!***routines called
!***end prologue       Input_Packed_Matrices_z
!
  SUBROUTINE Input_Packed_Matrices_z(file_directory,output_matrices)
  IMPLICIT NONE
  COMPLEX*16                            :: matrix_el
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: ij
  INTEGER                               :: ind_i
  INTEGER                               :: ind_j
  INTEGER                               :: index
  REAL*8                                :: gbytes
  INTEGER                               :: len
  INTEGER                               :: lenth
  CHARACTER(LEN=*)                      :: output_matrices
  CHARACTER(LEN=*)                      :: file_directory
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
      lwork = 10 * n3d
      write(iout,*) 'Diagonalizing Big Matrix'
      IF(non_orth) THEN
         ALLOCATE( h_mat_z(1:n3d,1:n3d), s_mat_z(1:n3d,1:n3d), eig(1:n3d),              &
                   work_z(1:lwork), rwork(1:lwork) ) 
         s_mat_z(1:n3d,1:n3d) = 0.d0
         Call Fill_Matrix_from_Disk(s_mat_z,non_zero_overlap_elements,50)
         h_mat_z(1:n3d,1:n3d) = 0.d0
         Call Fill_Matrix_from_Disk(h_mat_z,non_zero_hamiltonian_elements,50)
         call zhegv(l,'v','u',n3d,h_mat_z,n3d,s_mat_z,n3d,eig,work_z,                   &
                    lwork,rwork,info)
      ELSE
         ALLOCATE( h_mat_z(1:n3d,1:n3d), eig(1:n3d), work_z(1:lwork), rwork(1:lwork) )
         h_mat_z(1:n3d,1:n3d) = 0.d0
         Call Fill_Matrix_from_Disk(h_mat_z,non_zero_hamiltonian_elements,50)
         call zheev('v','u',n3d,h_mat_z,n3d,eig,work_z,lwork,rwork,info)
      END IF
      write(iout,*) 'Done Diagonalizing Big Matrix'
      local_title='eigenvalues'
      call prntcmn(local_title,eig,n3d,1,n3d,1,iout,'e')
      Call IOsys('write real eigenvalues to packed_matrices',n3d,eig,0,' ')
      Call IOsys('write real ground_state_eigenvector to packed_matrices',              &
                  2*n3d,h_mat_z,0,' ')
      Call IOsys('close packed_matrices',0,0,0,' ')
      Write(iout,*) 'Finished Diagonalizing Input Matrix: Quit'
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
        write(iout,*) 'Pack the Overlap Matrix'
        write(iout,*) 'Number Non_Zero Overlap In = ',non_zero_overlap_elements
        ALLOCATE(ibuf(2,non_zero_overlap_elements),                                     &
                 h_buf_z(non_zero_overlap_elements),diag_z(n3d))
        drop_tol=drop_overlap
        Call Fill_Buffers_from_Disk(ibuf,h_buf_z,diag_z,drop_tol,                       &
                                    non_zero_overlap_elements,50)
        write(iout,*) 'Number Non_Zero Overlap Out = ',non_zero_overlap_elements
        write(iout,*) 'Write the Overlap Matrix to Disk'
        Call Read_and_Write_Packed_Matrices(ibuf,h_buf_z,diag_z,                        &
                                            non_zero_overlap_elements,'overlap','write')
!       
!              Cholesky Decomposition and write to DIsk.
!
        write(iout,*)
        write(iout,*)'Form Cholesky Decomposition'
        ALLOCATE(upper_z(1:tri_size))
        Call Fill_Matrix_from_Buffers(upper_z,ibuf,h_buf_z,diag_z,                      &
                                      non_zero_overlap_elements)
        call zpptrf('u',upper_z,n3d,info)
        IF(info /= 0) THEN
           write(iout,*) 'Singular Overlap Matrix:Quit'
           Call lnkerr('Quit.  Singular Overlap Matrix')
        END IF
        write(iout,*) 'Opening File to Write Full Cholesky Decomposition'
        Call IOsys('open cholesky_decomposition as new',0,0,0,                          &
                    file_directory(1:len)//'/'//'Cholesky_Decomposition')
        write(iout,*)'Write Full Cholesky Decomposition to Disk'
        Call IOsys('write real cholesky_factor to cholesky_decomposition',              &
                    2*tri_size,upper_z,0,' ')
        write(iout,*) 'Closing File Holding Full Cholesky Decomposition'
        Call IOsys('close cholesky_decomposition',0,0,0,' ')
        drop_tol=drop_cholesky
        ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                      &
                 mat_var(1)%row_index(lenbuf),                                          &
                 mat_var(1)%packed_columns_z(lenbuf) ,                                  &
                 mat_var(1)%matrix_diagonal_z(n3d) )
        Call Pack_Triangle(upper_z,'upper')
        DEALLOCATE(mat_var(1)%non_zero_columns,                                         &
                   mat_var(1)%row_index,                                                &
                   mat_var(1)%packed_columns_z )
        ALLOCATE(lower_z(1:tri_size))
        Call Upper_to_Lower(upper_z,lower_z)
        ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                      &
                 mat_var(2)%row_index(lenbuf),                                          &
                 mat_var(2)%packed_columns_z(lenbuf) )
        Call Pack_Triangle(lower_z,'lower')
        DEALLOCATE(mat_var(2)%non_zero_columns,                                         &
                   mat_var(2)%row_index,                                                &
                   mat_var(2)%packed_columns_z,                                         &
                   mat_var(1)%matrix_diagonal_z )
        DEALLOCATE(ibuf, h_buf_z, lower_z)
        non_zero_cholesky_elements = tri_size
        ALLOCATE(ibuf(2,tri_size), h_buf_z(tri_size) )
        Call Fill_Buffers_from_Matrix(upper_z,ibuf,h_buf_z,diag_z,drop_tol,             &
                          non_zero_cholesky_elements)
        DEALLOCATE(upper_z)
        write(iout,*)
        write(iout,*) 'Pack the Cholesky Factors'
        write(iout,*) 'Number Non_Zero Cholesky Out = ', non_zero_cholesky_elements 
        write(iout,*)'Write Packed Cholesky Factors to Disk'
        Call Read_and_Write_Packed_Matrices(ibuf,h_buf_z,                               &
                                  diag_z,non_zero_cholesky_elements,                    &
                                  'cholesky','write')
        DEALLOCATE(ibuf, h_buf_z, diag_z )
     END IF
     READ(50) non_zero_hamiltonian_elements
     write(iout,*)
     write(iout,*) 'Pack the Hamiltonian Matrix'
     write(iout,*) 'Number Non_Zero Hamiltonian In = ',non_zero_hamiltonian_elements
     drop_tol=drop_hamiltonian
     Call Fill_Buffers_from_Disk(ibuf,h_buf_z,diag_z,drop_tol,                          &
                                 non_zero_hamiltonian_elements,50)
     write(iout,*) 'Number Non_Zero Hamiltonian Out = ',non_zero_hamiltonian_elements
     write(iout,*) 'Write Packed Hamiltonian to Disk'
     Call Read_and_Write_Packed_Matrices(ibuf,h_buf_z,diag_z,                           &
                                    non_zero_hamiltonian_elements,'hamiltonian','write')
     DEALLOCATE( ibuf, h_buf_z, diag_z )
  ELSE
!
!       Triangular Packing by Columns.
!
     IF(non_orth) THEN      
        ALLOCATE(triangle_overlap_z(tri_size))
        triangle_overlap_z(1:tri_size) = 0.d0
        Call Fill_Matrix_from_Disk(triangle_overlap_z,non_zero_overlap_elements,50)
        write(iout,*) 'Number Non_Zero Overlap In = ',non_zero_overlap_elements
        drop_tol=drop_overlap
        write(iout,*)
        write(iout,*) 'Pack the Overlap'
        ALLOCATE(mat_var(4)%non_zero_columns(n3d),                                      &
                 mat_var(4)%row_index(lenbuf),                                          &
                 mat_var(4)%packed_columns_z(lenbuf) ,                                  &
                 mat_var(4)%matrix_diagonal_z(n3d) )
        Call Pack_Triangle(triangle_overlap_z,'overlap')
        DEALLOCATE(mat_var(4)%non_zero_columns,                                         &
                   mat_var(4)%row_index,                                                &
                   mat_var(4)%packed_columns_z,                                         &
                   mat_var(4)%matrix_diagonal_z )
        ALLOCATE(upper_z(tri_size))
        upper_z(1:tri_size) = triangle_overlap_z(1:tri_size)
        write(iout,*)
        write(iout,*)'Form Cholesky Decomposition'
        call dpptrf('u',n3d,upper_z,info)
        IF(info /= 0) THEN
           write(iout,*) 'Singular Overlap Matrix:Quit'
           Call lnkerr('Quit.  Singular Overlap Matrix')
        END IF
        write(iout,*)
        write(iout,*) 'Pack the Cholesky Factors'
        drop_tol=drop_cholesky
        ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                      &
                 mat_var(1)%row_index(lenbuf),                                          &
                 mat_var(1)%packed_columns_z(lenbuf) ,                                  &
                 mat_var(1)%matrix_diagonal_z(n3d) )
        Call Pack_Triangle(upper_z,'upper')
        DEALLOCATE(mat_var(1)%non_zero_columns,                                         &
                   mat_var(1)%row_index,                                                &
                   mat_var(1)%packed_columns_z )
        Call Upper_to_Lower(upper_z,triangle_overlap_z)
        ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                      &
                 mat_var(2)%row_index(lenbuf),                                          &
                 mat_var(2)%packed_columns_z(lenbuf) )
        Call Pack_Triangle(triangle_overlap_z,'lower')
        DEALLOCATE(mat_var(2)%non_zero_columns,                                         &
                   mat_var(2)%row_index,                                                &
                   mat_var(2)%packed_columns_z,                                         &
                   mat_var(1)%matrix_diagonal_z )
        DEALLOCATE(triangle_overlap_z, upper_z)
     END IF
     ALLOCATE(triangle_hamiltonian_z(tri_size))
     triangle_hamiltonian_z(1:tri_size) = 0.d0
     Call Fill_Matrix_from_Disk(triangle_hamiltonian_z,non_zero_hamiltonian_elements,50)
     write(iout,*) 'Number Non_Zero Hamiltonian In = ',non_zero_hamiltonian_elements
     write(iout,*)
     write(iout,*) 'Pack the Hamiltonian'
     ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                         &
              mat_var(3)%row_index(lenbuf),                                             &
              mat_var(3)%packed_columns_z(lenbuf),                                      &
              mat_var(3)%matrix_diagonal_z(n3d) )
     drop_tol=drop_hamiltonian
     Call Pack_Triangle_z(triangle_hamiltonian_z,'hamiltonian')
     DEALLOCATE(mat_var(3)%non_zero_columns,                                            &
                mat_var(3)%row_index,                                                   &
                mat_var(3)%packed_columns_z,                                            &
                mat_var(3)%matrix_diagonal_z )
     DEALLOCATE(triangle_hamiltonian_z)
  END IF
!***********************************************************************
  END SUBROUTINE Input_Packed_Matrices_z
!***********************************************************************
!***********************************************************************
  END  MODULE Packed_Matrix_Module
!***********************************************************************
!***********************************************************************
