!***********************************************************************
                           MODULE Packed_Matrix_Module
                           USE Matrix_Print
                           USE Atomic_Matrices
                           USE Global_Time_Propagation_Module
                           USE Pack_Hamiltonian_Module
                           USE Preconditioner_Module
                           USE FEDVR_Shared, ONLY : file_loc
!                           USE input_output
!
                           IMPLICIT NONE
            CHARACTER(LEN=80)                    :: local_title
            CHARACTER(LEN=160)                   :: hold_str
            INTEGER, DIMENSION(5)                :: len
            TYPE(REAL_MATRIX)                    :: type_real_matrix
            TYPE(REAL_VECTOR)                    :: type_real_vector
            TYPE(COMPLEX_VECTOR)                 :: type_complex_vector
            TYPE(COMPLEX_MATRIX)                 :: type_complex_matrix
!***********************************************************************
!***********************************************************************
!
                          INTERFACE Packed_Matrix
                    MODULE PROCEDURE Packed_Matrix_d,                                    &
                                    Packed_Matrix_z
                          END INTERFACE Packed_Matrix
!
                          INTERFACE Input_Packed_IOsys_Matrix
                    MODULE PROCEDURE Input_Packed_IOsys_Matrix_d,                        &
                                    Input_Packed_IOsys_Matrix_z
                          END INTERFACE Input_Packed_IOsys_Matrix
!
                          INTERFACE Input_Packed_Matrices
                    MODULE PROCEDURE Input_Packed_Matrices_d,                            &
                                    Input_Packed_Matrices_z
                          END INTERFACE Input_Packed_Matrices
!
                          INTERFACE Read_Input_Matrices
                    MODULE PROCEDURE Read_Input_Matrices_d,                              &
                                     Read_Input_Matrices_z
                          END INTERFACE Read_Input_Matrices
!
                          INTERFACE Diagonalize_Packed_Matrices
                    MODULE PROCEDURE Diagonalize_Packed_Matrices_d,                      &
                                     Diagonalize_Packed_Matrices_z
                          END INTERFACE Diagonalize_Packed_Matrices
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                               Contains
!***********************************************************************
!***********************************************************************
!deck Packed_Matrix_d
!***begin prologue     Packed_Matrix_d
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
!***end prologue       Packed_Matrix_d
!
  SUBROUTINE Packed_Matrix_d(mat_d)
  IMPLICIT NONE
  TYPE(REAL_MATRICES)                   :: mat_d
  INTEGER                               :: lenth
  INTEGER                               :: IOSTAT
  REAL(isp)                             :: secnds
  REAL(isp)                             :: del_t
!
!
!*********************************************************************************************
!
!          Either reopen an already packed file or make a new packed file
!
!*********************************************************************************************
  IF (input_matrices(11:22) == 'iosys_format') THEN  ! Matrices have already been packed.  Just read them in.
      write(iout,*)
      write(iout,*) '                    *** Reopening File to Read Packed Matrices ***'
      write(iout,*)
      time(1) = secnds(0.0)
      file_loc=File_Directory(4)(1:len_dir(4))//'/'//'packed_matrix_file'
      Call pakstr(file_loc,len(1))
      Call IOsys('open packed_matrices as old',0,0,0,file_loc(1:len(1)))
      Call Input_Packed_IOsys_Matrix(mat_d)
      write(iout,*) '*** Closing Packed Matrix File ***'
      Call IOsys('close packed_matrices',0,0,0,' ')
      time(2) = secnds(0.0)
      del_t = time(2) - time(1)
      WRITE(iout,1)
      WRITE(iout,3) del_t
      WRITE(iout,1)
  ELSE                                               ! Pack the matrices in IOsys format 
      write(iout,*)
      write(iout,*) '          Reformatting Packed Matrices to IOsys Format'
      write(iout,*)
!
!        In this branch of the if, we come in with a packed matrix from the Drake codes, 
!        and go out with a packed matrix in IOsys format.  The packed matrix from IOsys 
!        is done in the standard or triangular manner.  In the standard manner an array of 
!        two indices and an array of elements that correspond to those 
!        indices is output along with the value of the number of non zero elements.  
!        Only the lower triangle is packed.  The other alternative is a packing of the 
!        triangle of the matrix.  This column oriented packing requires more extensive
!        logic but is needed to efficiently solve the linear equations when a 
!        non-orthogonal basis is used.
      time(1) = secnds(0.0) 
      write(iout,*)
      write(iout,*) '                    *** Opening File to Write Packed Matrices ***'
      write(iout,*)
      file_loc=File_Directory(4)(1:len_dir(4))//'/'//'packed_matrix_file'
      write(iout,*) file_loc
      Call pakstr(file_loc,len(1))
      Call IOsys('open packed_matrices as new',0,0,0,file_loc(1:len(1)))
      IF (mat%device == 'from_input' ) THEN
          Call Read_Input_Matrices(mat_d)
      ELSE
          Call pakstr(mat%device,len(1))    
          Call pakstr(species,len(2))    
          WRITE(iout,*) '        1. Hamiltonian File Name = '//mat%device(1:len(1))
          file_loc = File_Directory(4)(1:len_dir(4))//'/'//species(1:len(2))//'/'//mat%device(1:len(1))
          write(iout,*) file_loc
          Call pakstr(file_loc,len(3))
          write(iout,*) '        2. File Location = ',file_loc(1:len(3))
          write(iout,*)
          write(iout,*) '                    *** Processing Input Matrix File ***'
          write(iout,*)
          OPEN(UNIT=50,FILE=file_loc(1:len(3)),ACCESS='sequential', FORM='unformatted',   &
               IOSTAT=IOSTAT,STATUS='old')
          READ(50)
          Call Input_Packed_Matrices(mat_d)
          write(iout,*)
          write(iout,*) '                    *** Closing Input Matrix File ***'
          write(iout,*)
          CLOSE(50)
      END IF
      write(iout,*)
      write(iout,*) '                    *** Ending Writing of Packed Matrices ***'
      write(iout,*)
      Call IOsys ('close packed_matrices',0,0,0,' ')
      IF( reformat_input_matrix_only) THEN
!            We can opt to quit now if we wish or to re-read in the matrices for
!            a full calculation.
!
          write(iout,*) '                    *** Quitting After Reformat ***'
          time(2) = secnds(0.0)
          del_t = time(2) - time(1)
          WRITE(iout,1)
          WRITE(iout,2) del_t
          WRITE(iout,1)
          stop
      ELSE
          write(iout,*)
          write(iout,*) '                    *** Reopening File to Read Packed Matrices ***'
          write(iout,*)
          file_loc=File_Directory(4)(1:len_dir(4))//'/'//'packed_matrix_file'
          Call pakstr(file_loc,len(1))
          Call IOsys('open packed_matrices as old',0,0,0,file_loc(1:len(1))) 
          Call Input_Packed_IOsys_Matrix(mat_d)
          write(iout,*)
          write(iout,*) '               *** Closing File Packed Matrices ***'
          write(iout,*)
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
END SUBROUTINE Packed_Matrix_d
!***********************************************************************
!***********************************************************************
!deck Packed_Matrix_z
!***begin prologue     Packed_Matrix_z
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
!***end prologue       Packed_Matrix_z
!
  SUBROUTINE Packed_Matrix_z(mat_z)
  IMPLICIT NONE
  TYPE(COMPLEX_MATRICES)                :: mat_z
  INTEGER                               :: lenth
  INTEGER                               :: IOSTAT
  REAL(isp)                             :: secnds
  REAL(isp)                             :: del_t
!
!
!*********************************************************************************************
!
!          Either reopen an already packed file or make a new packed file
!
!*********************************************************************************************
  IF (input_matrices(11:22) == 'iosys_format') THEN  ! Matrices have already been packed.  Just read them in.
      write(iout,*)
      write(iout,*) '                    *** Reopening File to Read Packed Matrices ***'
      write(iout,*)
      time(1) = secnds(0.0)
      file_loc=File_Directory(4)(1:len_dir(4))//'/'//'packed_matrix_file'
      Call pakstr(file_loc,len(1))
      Call IOsys('open packed_matrices as old',0,0,0,file_loc(1:len(1)))
      Call Input_Packed_IOsys_Matrix(mat_z)
      write(iout,*) '*** Closing Packed Matrix File ***'
      Call IOsys('close packed_matrices',0,0,0,' ')
      time(2) = secnds(0.0)
      del_t = time(2) - time(1)
      WRITE(iout,1)
      WRITE(iout,3) del_t
      WRITE(iout,1)
  ELSE                                               ! Pack the matrices in IOsys format 
      write(iout,*)
      write(iout,*) '          Reformatting Packed Matrices to IOsys Format'
      write(iout,*)
!
!        In this branch of the if, we come in with a packed matrix from the Drake codes, 
!        and go out with a packed matrix in IOsys format.  The packed matrix from IOsys 
!        is done in the standard or triangular manner.  In the standard manner an array of 
!        two indices and an array of elements that correspond to those 
!        indices is output along with the value of the number of non zero elements.  
!        Only the lower triangle is packed.  The other alternative is a packing of the 
!        triangle of the matrix.  This column oriented packing requires more extensive
!        logic but is needed to efficiently solve the linear equations when a 
!        non-orthogonal basis is used.
      time(1) = secnds(0.0) 
      write(iout,*)
      write(iout,*) '                    *** Opening File to Write Packed Matrices ***'
      write(iout,*)
      file_loc=File_Directory(4)(1:len_dir(4))//'/'//'packed_matrix_file'
      write(iout,*) file_loc
      Call pakstr(file_loc,len(1))
      Call IOsys('open packed_matrices as new',0,0,0,file_loc(1:len(1)))
      IF (mat%device == 'from_input' ) THEN
          Call Read_Input_Matrices(mat_z)
      ELSE
          Call pakstr(mat%device,len(1))    
          Call pakstr(species,len(2))    
          WRITE(iout,*) '        1. Hamiltonian File Name = '//mat%device(1:len(1))
          file_loc = File_Directory(4)(1:len_dir(4))//'/'//species(1:len(2))//'/'//mat%device(1:len(1))
          write(iout,*) file_loc
          Call pakstr(file_loc,len(3))
          write(iout,*) '        2. File Location = ',file_loc(1:len(3))
          write(iout,*)
          write(iout,*) '                    *** Processing Input Matrix File ***'
          write(iout,*)
          OPEN(UNIT=50,FILE=file_loc(1:len(3)),ACCESS='sequential', FORM='unformatted',   &
               IOSTAT=IOSTAT,STATUS='old')
          READ(50)
          Call Input_Packed_Matrices(mat_z)
          write(iout,*)
          write(iout,*) '                    *** Closing Input Matrix File ***'
          write(iout,*)
          CLOSE(50)
      END IF
      write(iout,*)
      write(iout,*) '                    *** Ending Writing of Packed Matrices ***'
      write(iout,*)
      Call IOsys ('close packed_matrices',0,0,0,' ')
      IF( reformat_input_matrix_only) THEN
!            We can opt to quit now if we wish or to re-read in the matrices for
!            a full calculation.
!
          write(iout,*) '                    *** Quitting After Reformat ***'
          time(2) = secnds(0.0)
          del_t = time(2) - time(1)
          WRITE(iout,1)
          WRITE(iout,2) del_t
          WRITE(iout,1)
          stop
      ELSE
          write(iout,*)
          write(iout,*) '                    *** Reopening File to Read Packed Matrices ***'
          write(iout,*)
          file_loc=File_Directory(4)(1:len_dir(4))//'/'//'packed_matrix_file'
          Call pakstr(file_loc,len(1))
          Call IOsys('open packed_matrices as old',0,0,0,file_loc(1:len(1))) 
          Call Input_Packed_IOsys_Matrix(mat_z)
          write(iout,*)
          write(iout,*) '               *** Closing File Packed Matrices ***'
          write(iout,*)
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
END SUBROUTINE Packed_Matrix_z
!***********************************************************************
!***********************************************************************
!deck Input_Packed_IOsys_Matrix_d
!***begin prologue     Input_Packed_IOsys_Matrix_d
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
!***end prologue       Input_Packed_IOsys_Matrix_d
!
  SUBROUTINE Input_Packed_IOsys_Matrix_d(mat_d)
  IMPLICIT NONE
  TYPE(REAL_MATRICES)                   :: mat_d
  INTEGER                               :: max_buf
!
!
  mat%non_zero_overlap_elements = 0
  mat%non_zero_cholesky_elements = 0
  write(iout,*)
  write(iout,*) '                    *** Reading in Packed Matrices in Triangle Form ***'
  write(iout,*)
!
!    The matrices here are already packed in column form and only need to be read in
!
  IF (non_orth) THEN
      IF (overlap_to_disk) THEN   
          write(iout,*)
          write(iout,*) '                    *** Reading Packed Overlap ***'
          write(iout,*)
          Call IOsys('read integer "number_non_zero_overlap_elements" from '//      &
                     'packed_matrices',1,mat_var(1)%number,0,' ')
          ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                &
                   mat_var(1)%row_index(mat_var(1)%number),                         &
                   mat_var(1)%mat_d%packed_columns(mat_var(1)%number),              &
                   mat_var(1)%mat_d%matrix_diagonal(n3d) )
          Call Read_and_Write_Column_Packed_Matrices                                &
                                           (mat_var(1)%mat_d%packed_columns,        &
                                            mat_var(1)%non_zero_columns,            &
                                            mat_var(1)%row_index,                   &
                                            'read',                                 &
                                            'overlap',                              &
                                             mat_var(1)%number,                     &
                                             mat_var(1)%mat_d%matrix_diagonal )
      END IF
      write(iout,*)
      write(iout,*) '                    *** Reading Packed Upper Cholesky Factor ***'
      write(iout,*)
      Call IOsys('read integer "number_non_zero_upper_elements" from '//        &
                 'packed_matrices',1,mat_var(1)%number,0,' ')
      ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                &
               mat_var(1)%row_index(mat_var(1)%number),                         &
               mat_var(1)%mat_d%packed_columns(mat_var(1)%number),              &
               mat_var(1)%mat_d%matrix_diagonal(n3d) )
      Call Read_and_Write_Column_Packed_Matrices                                &
                                       (mat_var(1)%mat_d%packed_columns,        &
                                        mat_var(1)%non_zero_columns,            &
                                        mat_var(1)%row_index,                   &
                                        'read',                                 &
                                        'upper',                                &
                                        mat_var(1)%number,                      &
                                        mat_var(1)%mat_d%matrix_diagonal )
      write(iout,*) '          *** Reading Packed Lower Cholesky Factor ***'
      Call IOsys('read integer "number_non_zero_lower_elements" from '//        &
                 'packed_matrices',1,mat_var(2)%number,0,' ')
      ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                &
               mat_var(2)%row_index(mat_var(2)%number),                         &
               mat_var(2)%mat_d%packed_columns(mat_var(2)%number) )
      Call Read_and_Write_Column_Packed_Matrices                                &
                                       (mat_var(2)%mat_d%packed_columns,        &
                                        mat_var(2)%non_zero_columns,            &
                                        mat_var(2)%row_index,                   &
                                        'read',                                 &
                                        'lower',                                &
                                        mat_var(2)%number )
  END IF
  write(iout,*)
  write(iout,*) '                    *** Reading Packed Hamiltonian ***'
  write(iout,*)
  Call IOsys('read integer "number_non_zero_hamiltonian_elements" from '//      &
             'packed_matrices',1,mat_var(3)%number,0,' ')
  ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                    &
           mat_var(3)%row_index(mat_var(3)%number),                             &
           mat_var(3)%mat_d%packed_columns(mat_var(3)%number),                  &
           mat_var(3)%mat_d%matrix_diagonal(n3d) )
  Call Read_and_Write_Column_Packed_Matrices                                    &
                                   (mat_var(3)%mat_d%packed_columns,            &
                                    mat_var(3)%non_zero_columns,                &
                                    mat_var(3)%row_index,                       &
                                   'read',                                      &
                                   'hamiltonian',                               &
                                    mat_var(3)%number,                          &
                                    mat_var(3)%mat_d%matrix_diagonal )
!***********************************************************************
  END SUBROUTINE Input_Packed_IOsys_Matrix_d
!***********************************************************************
!***********************************************************************
!deck Input_Packed_IOsys_Matrix_z
!***begin prologue     Input_Packed_IOsys_Matrix_z
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
!***end prologue       Input_Packed_IOsys_Matrix_z
!
  SUBROUTINE Input_Packed_IOsys_Matrix_z(mat_z)
  IMPLICIT NONE
  TYPE(COMPLEX_MATRICES)                :: mat_z
  INTEGER                               :: max_buf
!
!
  mat%non_zero_overlap_elements = 0
  mat%non_zero_cholesky_elements = 0
  write(iout,*)
  write(iout,*) '                    *** Reading in Packed Matrices in Triangle Form ***'
  write(iout,*)
!
!    The matrices here are already packed in column form and only need to be read in
!
  IF (non_orth) THEN
      IF (overlap_to_disk) THEN   
          write(iout,*)
          write(iout,*) '                    *** Reading Packed Overlap ***'
          write(iout,*)
          Call IOsys('read integer "number_non_zero_overlap_elements" from '//      &
                     'packed_matrices',1,mat_var(1)%number,0,' ')
          ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                &
                   mat_var(1)%row_index(mat_var(1)%number),                         &
                   mat_var(1)%mat_z%packed_columns(mat_var(1)%number),              &
                   mat_var(1)%mat_z%matrix_diagonal(n3d) )
          Call Read_and_Write_Column_Packed_Matrices                                &
                                           (mat_var(1)%mat_z%packed_columns,        &
                                            mat_var(1)%non_zero_columns,            &
                                            mat_var(1)%row_index,                   &
                                            'read',                                 &
                                            'overlap',                              &
                                             mat_var(1)%number,                     &
                                             mat_var(1)%mat_z%matrix_diagonal )
      END IF
      write(iout,*)
      write(iout,*) '                    *** Reading Packed Upper Cholesky Factor ***'
      write(iout,*)
      Call IOsys('read integer "number_non_zero_upper_elements" from '//        &
                 'packed_matrices',1,mat_var(1)%number,0,' ')
      ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                &
               mat_var(1)%row_index(mat_var(1)%number),                         &
               mat_var(1)%mat_z%packed_columns(mat_var(1)%number),              &
               mat_var(1)%mat_z%matrix_diagonal(n3d) )
      Call Read_and_Write_Column_Packed_Matrices                                &
                                       (mat_var(1)%mat_z%packed_columns,        &
                                        mat_var(1)%non_zero_columns,            &
                                        mat_var(1)%row_index,                   &
                                        'read',                                 &
                                        'upper',                                &
                                        mat_var(1)%number,                      &
                                        mat_var(1)%mat_z%matrix_diagonal )
      write(iout,*) '          *** Reading Packed Lower Cholesky Factor ***'
      Call IOsys('read integer "number_non_zero_lower_elements" from '//        &
                 'packed_matrices',1,mat_var(2)%number,0,' ')
      ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                &
               mat_var(2)%row_index(mat_var(2)%number),                         &
               mat_var(2)%mat_z%packed_columns(mat_var(2)%number) )
      Call Read_and_Write_Column_Packed_Matrices                                &
                                       (mat_var(2)%mat_z%packed_columns,        &
                                        mat_var(2)%non_zero_columns,            &
                                        mat_var(2)%row_index,                   &
                                        'read',                                 &
                                        'lower',                                &
                                        mat_var(2)%number )
  END IF
  write(iout,*)
  write(iout,*) '                    *** Reading Packed Hamiltonian ***'
  write(iout,*)
  Call IOsys('read integer "number_non_zero_hamiltonian_elements" from '//      &
             'packed_matrices',1,mat_var(3)%number,0,' ')
  ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                    &
           mat_var(3)%row_index(mat_var(3)%number),                             &
           mat_var(3)%mat_z%packed_columns(mat_var(3)%number),                  &
           mat_var(3)%mat_z%matrix_diagonal(n3d) )
  Call Read_and_Write_Column_Packed_Matrices                                    &
                                   (mat_var(3)%mat_z%packed_columns,            &
                                    mat_var(3)%non_zero_columns,                &
                                    mat_var(3)%row_index,                       &
                                   'read',                                      &
                                   'hamiltonian',                               &
                                    mat_var(3)%number,                          &
                                    mat_var(3)%mat_z%matrix_diagonal )
!***********************************************************************
  END SUBROUTINE Input_Packed_IOsys_Matrix_z
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
  SUBROUTINE Input_Packed_Matrices_d(mat_d)
  IMPLICIT NONE
  TYPE(REAL_MATRICES)                   :: mat_d
  REAL(idp)                             :: gbytes
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
  IF (diagonalize_only) THEN
      Call Diagonalize_Packed_Matrices  (mat_d,                                       &
                                        'disk',                                       &
                                         n3d,                                          &
                                         .false.)
      stop
  END IF
!
! This is the option to reformat the input matrix in IOsys format and then return 
! to the calling routine.

!
!    Triangular Packing by Columns.
!
  IF(non_orth) THEN      
     ALLOCATE(mat_d%triangle_overlap(tri_size))
     Call Fill_Matrix_from_Disk (                                                       &
                                 mat_d%triangle_overlap,                                &
                                 mat%non_zero_overlap_elements,                         &
                                 50)
     write(iout,*) '           Number Non_Zero Overlap In = ',mat%non_zero_overlap_elements
     IF (overlap_to_disk) THEN   
         drop_tol=drop_overlap
         write(iout,*)
         write(iout,*) '                    *** Pack the Overlap ***'
         write(iout,*)
         ALLOCATE(mat_var(4)%non_zero_columns(n3d),                                     &
                  mat_var(4)%row_index(lenbuf),                                         &
                  mat_var(4)%mat_d%packed_columns(lenbuf) ,                             &
                  mat_var(4)%mat_d%matrix_diagonal(n3d) )
         Call Pack_Triangle ('overlap',                                                 &
                              mat_d%triangle_overlap,                                   &
                              mat_var(4)%mat_d%packed_columns,                          &
                              mat_var(4)%non_zero_columns,                              &
                              mat_var(4)%row_index,                                     &
                              mat_var(4)%number,                                        &
                              mat_var(4)%mat_d%matrix_diagonal)
         DEALLOCATE(mat_var(4)%non_zero_columns,                                        &
                    mat_var(4)%row_index,                                               &
                    mat_var(4)%mat_d%packed_columns,                                    &
                    mat_var(4)%mat_d%matrix_diagonal )
     END IF
     ALLOCATE(mat_d%upper(tri_size))
     mat_d%upper(1:tri_size) = mat_d%triangle_overlap(1:tri_size)
     write(iout,*)
     write(iout,*)'                    *** Form Cholesky Decomposition ***'
     write(iout,*)
     call dpptrf('u',n3d,mat_d%upper,info)
     IF(info /= 0) THEN
        write(iout,*)
        write(iout,*) '                    $$$ Catastrophe - Singular Overlap Matrix:Quit $$$'
        write(iout,*)
        Call lnkerr('Quit.  Singular Overlap Matrix')
     END IF
     write(iout,*)
     write(iout,*) '                    *** Pack the Cholesky Factors ***'
     write(iout,*)
     drop_tol=drop_cholesky
     ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                         &
              mat_var(1)%row_index(lenbuf),                                             &
              mat_var(1)%mat_d%packed_columns(lenbuf) ,                                 &
              mat_var(1)%mat_d%matrix_diagonal(n3d) )
     Call Pack_Triangle ('upper',                                                       &
                          mat_d%upper,                                                  &
                          mat_var(1)%mat_d%packed_columns,                              &
                          mat_var(1)%non_zero_columns,                                  &
                          mat_var(1)%row_index,                                         &
                          mat_var(1)%number,                                            &
                          mat_var(1)%mat_d%matrix_diagonal)
     DEALLOCATE(mat_var(1)%non_zero_columns,                                            &
                mat_var(1)%row_index,                                                   &
                mat_var(1)%mat_d%packed_columns )
     Call Upper_to_Lower(mat_d%upper,mat_d%triangle_overlap)
     DEALLOCATE(mat_d%upper)
     ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                         &
              mat_var(2)%row_index(lenbuf),                                             &
              mat_var(2)%mat_d%packed_columns(lenbuf) )
     Call Pack_Triangle ('lower',                                                       &
                          mat_d%triangle_overlap,                                       &
                          mat_var(2)%mat_d%packed_columns,                              &
                          mat_var(2)%non_zero_columns,                                  &
                          mat_var(2)%row_index,                                         &
                          mat_var(2)%number,                                            &
                          mat_var(1)%mat_d%matrix_diagonal)
     DEALLOCATE(mat_var(2)%non_zero_columns,                                            &
                mat_var(2)%row_index,                                                   &
                mat_var(2)%mat_d%packed_columns,                                        &
                mat_var(1)%mat_d%matrix_diagonal )
     DEALLOCATE(mat_d%triangle_overlap)
  END IF
  ALLOCATE(mat_d%triangle_hamiltonian(1:tri_size))   
  Call Fill_Matrix_from_Disk (                                                          &
                              mat_d%triangle_hamiltonian,                               &
                              non_zero_hamiltonian_elements,                            &
                              50)
  write(iout,*)
  write(iout,*) '          Number Non_Zero Hamiltonian In = ',non_zero_hamiltonian_elements
  write(iout,*)
  write(iout,*) '                    *** Pack the Hamiltonian ***'
  ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                            &
           mat_var(3)%row_index(lenbuf),                                                &
           mat_var(3)%mat_d%packed_columns(lenbuf),                                     &
           mat_var(3)%mat_d%matrix_diagonal(n3d) )
  drop_tol=drop_hamiltonian
  Call Pack_Triangle ('hamiltonian',                                                    &
                       mat_d%triangle_hamiltonian,                                      &
                       mat_var(3)%mat_d%packed_columns,                                 &
                       mat_var(3)%non_zero_columns,                                     &
                       mat_var(3)%row_index,                                            &
                       mat_var(3)%number,                                               &
                       mat_var(3)%mat_d%matrix_diagonal)
  DEALLOCATE(mat_var(3)%non_zero_columns,                                               &
             mat_var(3)%row_index,                                                      &
             mat_var(3)%mat_d%packed_columns,                                           &
             mat_var(3)%mat_d%matrix_diagonal )
  DEALLOCATE(mat_d%triangle_hamiltonian)
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
  SUBROUTINE Input_Packed_Matrices_z(mat_z)
  IMPLICIT NONE
  TYPE(COMPLEX_MATRICES)                :: mat_z
  REAL(idp)                             :: gbytes
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
  IF (diagonalize_only) THEN
      Call Diagonalize_Packed_Matrices  (mat_z,                                        &
                                         'disk',                                       &
                                         n3d,                                          &
                                         .false.)
      stop
  END IF
!
! This is the option to reformat the input matrix in IOsys format and then return 
! to the calling routine.

!
!    Triangular Packing by Columns.
!
  IF(non_orth) THEN      
     ALLOCATE(mat_z%triangle_overlap(tri_size))
     Call Fill_Matrix_from_Disk (                                                       &
                                 mat_z%triangle_overlap,                                &
                                 mat%non_zero_overlap_elements,                         &
                                 50)
     write(iout,*) '           Number Non_Zero Overlap In = ',mat%non_zero_overlap_elements
     IF (overlap_to_disk) THEN   
         drop_tol=drop_overlap
         write(iout,*)
         write(iout,*) '                    *** Pack the Overlap ***'
         write(iout,*)
         ALLOCATE(mat_var(4)%non_zero_columns(n3d),                                     &
                  mat_var(4)%row_index(lenbuf),                                         &
                  mat_var(4)%mat_z%packed_columns(lenbuf) ,                             &
                  mat_var(4)%mat_z%matrix_diagonal(n3d) )
         Call Pack_Triangle ('overlap',                                                 &
                              mat_z%triangle_overlap,                                   &
                              mat_var(4)%mat_z%packed_columns,                          &
                              mat_var(4)%non_zero_columns,                              &
                              mat_var(4)%row_index,                                     &
                              mat_var(4)%number,                                        &
                              mat_var(4)%mat_z%matrix_diagonal)
         DEALLOCATE(mat_var(4)%non_zero_columns,                                        &
                    mat_var(4)%row_index,                                               &
                    mat_var(4)%mat_z%packed_columns,                                    &
                    mat_var(4)%mat_z%matrix_diagonal )
     END IF
     ALLOCATE(mat_z%upper(tri_size))
     mat_z%upper(1:tri_size) = mat_z%triangle_overlap(1:tri_size)
     write(iout,*)
     write(iout,*)'                    *** Form Cholesky Decomposition ***'
     write(iout,*)
     call zpptrf('u',n3d,mat_z%upper,info)
     IF(info /= 0) THEN
        write(iout,*)
        write(iout,*) '                    $$$ Catastrophe - Singular Overlap Matrix:Quit $$$'
        write(iout,*)
        Call lnkerr('Quit.  Singular Overlap Matrix')
     END IF
     write(iout,*)
     write(iout,*) '                    *** Pack the Cholesky Factors ***'
     write(iout,*)
     drop_tol=drop_cholesky
     ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                         &
              mat_var(1)%row_index(lenbuf),                                             &
              mat_var(1)%mat_z%packed_columns(lenbuf) ,                                 &
              mat_var(1)%mat_z%matrix_diagonal(n3d) )
     Call Pack_Triangle ('upper',                                                       &
                          mat_z%upper,                                                  &
                          mat_var(1)%mat_z%packed_columns,                              &
                          mat_var(1)%non_zero_columns,                                  &
                          mat_var(1)%row_index,                                         &
                          mat_var(1)%number,                                            &
                          mat_var(1)%mat_z%matrix_diagonal)
     DEALLOCATE(mat_var(1)%non_zero_columns,                                            &
                mat_var(1)%row_index,                                                   &
                mat_var(1)%mat_z%packed_columns )
     Call Upper_to_Lower(mat_z%upper,mat_z%triangle_overlap)
     DEALLOCATE(mat_z%upper)
     ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                         &
              mat_var(2)%row_index(lenbuf),                                             &
              mat_var(2)%mat_z%packed_columns(lenbuf) )
     Call Pack_Triangle ('lower',                                                       &
                          mat_z%triangle_overlap,                                       &
                          mat_var(2)%mat_z%packed_columns,                              &
                          mat_var(2)%non_zero_columns,                                  &
                          mat_var(2)%row_index,                                         &
                          mat_var(2)%number,                                            &
                          mat_var(1)%mat_z%matrix_diagonal)
     DEALLOCATE(mat_var(2)%non_zero_columns,                                            &
                mat_var(2)%row_index,                                                   &
                mat_var(2)%mat_z%packed_columns,                                        &
                mat_var(1)%mat_z%matrix_diagonal )
     DEALLOCATE(mat_z%triangle_overlap)
  END IF
  ALLOCATE(mat_z%triangle_hamiltonian(1:tri_size))   
  Call Fill_Matrix_from_Disk (                                                          &
                              mat_z%triangle_hamiltonian,                               &
                              mat%non_zero_hamiltonian_elements,                        &
                              50)
  write(iout,*)
  write(iout,*) '          Number Non_Zero Hamiltonian In = ',                          &
                           mat%non_zero_hamiltonian_elements
  write(iout,*)
  write(iout,*) '                    *** Pack the Hamiltonian ***'
  ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                            &
           mat_var(3)%row_index(lenbuf),                                                &
           mat_var(3)%mat_z%packed_columns(lenbuf),                                     &
           mat_var(3)%mat_z%matrix_diagonal(n3d) )
  drop_tol=drop_hamiltonian
  Call Pack_Triangle ('hamiltonian',                                                    &
                       mat_z%triangle_hamiltonian,                                      &
                       mat_var(3)%mat_z%packed_columns,                                 &
                       mat_var(3)%non_zero_columns,                                     &
                       mat_var(3)%row_index,                                            &
                       mat_var(3)%number,                                               &
                       mat_var(3)%mat_z%matrix_diagonal)
  DEALLOCATE(mat_var(3)%non_zero_columns,                                               &
             mat_var(3)%row_index,                                                      &
             mat_var(3)%mat_z%packed_columns,                                           &
             mat_var(3)%mat_z%matrix_diagonal )
  DEALLOCATE(mat_z%triangle_hamiltonian)
!***********************************************************************
  END SUBROUTINE Input_Packed_Matrices_z
!***********************************************************************
!***********************************************************************
!Deck Read_Input_Matrices_d
!***begin prologue     Read_Input_Matrices_d
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
!***end prologue       Read_Input_Matrices_d
!
  SUBROUTINE Read_Input_Matrices_d(mat_d)
  IMPLICIT NONE
  TYPE(REAL_MATRICES)                     :: mat_d
!
! This is the option to reformat the input matrix in IOsys format and then return 
! to the calling routine.

!
!    Triangular Packing by Columns.
!
  IF(non_orth) THEN      
     ALLOCATE(mat_d%triangle_overlap(tri_size))
     Call Fill_Matrix_from_Input ('$triangle_overlap',                              &
                                 mat_d%triangle_overlap,n3d)
     drop_tol=drop_overlap
     write(iout,*)
     write(iout,*) '                    *** Pack the Overlap ***'
     write(iout,*)
     ALLOCATE(mat_var(4)%non_zero_columns(n3d),                                     &
              mat_var(4)%row_index(lenbuf),                                         &
              mat_var(4)%mat_d%packed_columns(lenbuf) ,                             &
              mat_var(4)%mat_d%matrix_diagonal(n3d) )
     Call Pack_Triangle ('overlap',                                                 &
                          mat_d%triangle_overlap,                                   &
                          mat_var(4)%mat_d%packed_columns,                          &
                          mat_var(4)%non_zero_columns,                              &
                          mat_var(4)%row_index,                                     &
                          mat_var(4)%number,                                        &
                          mat_var(4)%mat_d%matrix_diagonal)
     DEALLOCATE(mat_var(4)%non_zero_columns,                                        &
                mat_var(4)%row_index,                                               &
                mat_var(4)%mat_d%packed_columns,                                    &
                mat_var(4)%mat_d%matrix_diagonal )
     ALLOCATE(mat_d%upper(tri_size))
     mat_d%upper(1:tri_size) = mat_d%triangle_overlap(1:tri_size)
     write(iout,*)
     write(iout,*)'                    ***Form Cholesky Decomposition ***'
     write(iout,*)
     call dpptrf('u',n3d,mat_d%upper,info)
     IF(info /= 0) THEN
        write(iout,*)
        write(iout,*) '                    $$$ Catastrophe - Singular Overlap Matrix:Quit $$$'
        write(iout,*)
        Call lnkerr('Quit.  Singular Overlap Matrix')
     END IF
     write(iout,*)
     write(iout,*) '                    *** Pack the Cholesky Factors ***'
     write(iout,*)
     drop_tol=drop_cholesky
     ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                         &
              mat_var(1)%row_index(lenbuf),                                             &
              mat_var(1)%mat_d%packed_columns(lenbuf) ,                                 &
              mat_var(1)%mat_d%matrix_diagonal(n3d) )
     Call Pack_Triangle ('upper',                                                       &
                          mat_d%upper,                                                  &
                          mat_var(1)%mat_d%packed_columns,                              &
                          mat_var(1)%non_zero_columns,                                  &
                          mat_var(1)%row_index,                                         &
                          mat_var(1)%number,                                            &
                          mat_var(1)%mat_d%matrix_diagonal)
     DEALLOCATE(mat_var(1)%non_zero_columns,                                            &
                mat_var(1)%row_index,                                                   &
                mat_var(1)%mat_d%packed_columns )
     Call Upper_to_Lower(mat_d%upper,mat_d%triangle_overlap)
     DEALLOCATE(mat_d%upper)
     ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                         &
              mat_var(2)%row_index(lenbuf),                                             &
              mat_var(2)%mat_d%packed_columns(lenbuf) )
     Call Pack_Triangle ('lower',                                                       &
                          mat_d%triangle_overlap,                                       &
                          mat_var(2)%mat_d%packed_columns,                              &
                          mat_var(2)%non_zero_columns,                                  &
                          mat_var(2)%row_index,                                         &
                          mat_var(2)%number,                                            &
                          mat_var(1)%mat_d%matrix_diagonal)
     DEALLOCATE(mat_var(2)%non_zero_columns,                                            &
                mat_var(2)%row_index,                                                   &
                mat_var(2)%mat_d%packed_columns,                                        &
                mat_var(1)%mat_d%matrix_diagonal )
  END IF
  ALLOCATE(mat_d%triangle_hamiltonian(1:tri_size))
  Call Fill_Matrix_from_Input ('$triangle_hamiltonian',                                 &
                              mat_d%triangle_hamiltonian,n3d)
  write(iout,*)
  write(iout,*) '          *** Pack the Hamiltonian ***'
  ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                            &
           mat_var(3)%row_index(lenbuf),                                                &
           mat_var(3)%mat_d%packed_columns(lenbuf),                                     &
           mat_var(3)%mat_d%matrix_diagonal(n3d) )
  drop_tol=drop_hamiltonian
  Call Pack_Triangle ('hamiltonian',                                                    &
                       mat_d%triangle_hamiltonian,                                      &
                       mat_var(3)%mat_d%packed_columns,                                 &
                       mat_var(3)%non_zero_columns,                                     &
                       mat_var(3)%row_index,                                            &
                       mat_var(3)%number,                                               &
                       mat_var(3)%mat_d%matrix_diagonal)
  DEALLOCATE(mat_var(3)%non_zero_columns,                                               &
             mat_var(3)%row_index,                                                      &
             mat_var(3)%mat_d%packed_columns,                                           &
             mat_var(3)%mat_d%matrix_diagonal )
!***********************************************************************
  END SUBROUTINE Read_Input_Matrices_d
!***********************************************************************
!***********************************************************************
!Deck Read_Input_Matrices_z
!***begin prologue     Read_Input_Matrices_z
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
!***end prologue       Read_Input_Matrices_z
!
  SUBROUTINE Read_Input_Matrices_z(mat_z)
  IMPLICIT NONE
  TYPE(COMPLEX_MATRICES)                     :: mat_z
!
! This is the option to reformat the input matrix in IOsys format and then return 
! to the calling routine.

!
!    Triangular Packing by Columns.
!
  IF(non_orth) THEN      
     ALLOCATE(mat_z%triangle_overlap(tri_size))
     Call Fill_Matrix_from_Input ('$triangle_overlap',                              &
                                 mat_z%triangle_overlap,n3d)
     drop_tol=drop_overlap
     write(iout,*)
     write(iout,*) '                    *** Pack the Overlap ***'
     write(iout,*)
     ALLOCATE(mat_var(4)%non_zero_columns(n3d),                                     &
              mat_var(4)%row_index(lenbuf),                                         &
              mat_var(4)%mat_z%packed_columns(lenbuf) ,                             &
              mat_var(4)%mat_z%matrix_diagonal(n3d) )
     Call Pack_Triangle ('overlap',                                                 &
                          mat_z%triangle_overlap,                                   &
                          mat_var(4)%mat_z%packed_columns,                          &
                          mat_var(4)%non_zero_columns,                              &
                          mat_var(4)%row_index,                                     &
                          mat_var(4)%number,                                        &
                          mat_var(4)%mat_z%matrix_diagonal)
     DEALLOCATE(mat_var(4)%non_zero_columns,                                        &
                mat_var(4)%row_index,                                               &
                mat_var(4)%mat_z%packed_columns,                                    &
                mat_var(4)%mat_z%matrix_diagonal )
     ALLOCATE(mat_z%upper(tri_size))
     mat_z%upper(1:tri_size) = mat_z%triangle_overlap(1:tri_size)
     write(iout,*)
     write(iout,*)'                    ***Form Cholesky Decomposition ***'
     write(iout,*)
     call zpptrf('u',n3d,mat_z%upper,info)
     IF(info /= 0) THEN
        write(iout,*)
        write(iout,*) '                    $$$ Catastrophe - Singular Overlap Matrix:Quit $$$'
        write(iout,*)
        Call lnkerr('Quit.  Singular Overlap Matrix')
     END IF
     write(iout,*)
     write(iout,*) '                    *** Pack the Cholesky Factors ***'
     write(iout,*)
     drop_tol=drop_cholesky
     ALLOCATE(mat_var(1)%non_zero_columns(n3d),                                         &
              mat_var(1)%row_index(lenbuf),                                             &
              mat_var(1)%mat_z%packed_columns(lenbuf) ,                                 &
              mat_var(1)%mat_z%matrix_diagonal(n3d) )
     Call Pack_Triangle ('upper',                                                       &
                          mat_z%upper,                                                  &
                          mat_var(1)%mat_z%packed_columns,                              &
                          mat_var(1)%non_zero_columns,                                  &
                          mat_var(1)%row_index,                                         &
                          mat_var(1)%number,                                            &
                          mat_var(1)%mat_z%matrix_diagonal)
     DEALLOCATE(mat_var(1)%non_zero_columns,                                            &
                mat_var(1)%row_index,                                                   &
                mat_var(1)%mat_z%packed_columns )
     Call Upper_to_Lower(mat_z%upper,mat_z%triangle_overlap)
     DEALLOCATE(mat_z%upper)
     ALLOCATE(mat_var(2)%non_zero_columns(n3d),                                         &
              mat_var(2)%row_index(lenbuf),                                             &
              mat_var(2)%mat_z%packed_columns(lenbuf) )
     Call Pack_Triangle ('lower',                                                       &
                          mat_z%triangle_overlap,                                       &
                          mat_var(2)%mat_z%packed_columns,                              &
                          mat_var(2)%non_zero_columns,                                  &
                          mat_var(2)%row_index,                                         &
                          mat_var(2)%number,                                            &
                          mat_var(1)%mat_z%matrix_diagonal)
     DEALLOCATE(mat_var(2)%non_zero_columns,                                            &
                mat_var(2)%row_index,                                                   &
                mat_var(2)%mat_z%packed_columns,                                        &
                mat_var(1)%mat_z%matrix_diagonal )
  END IF
  ALLOCATE(mat_z%triangle_hamiltonian(1:tri_size))
  Call Fill_Matrix_from_Input ('$triangle_hamiltonian',                                 &
                              mat_z%triangle_hamiltonian,n3d)
  write(iout,*)
  write(iout,*) '          *** Pack the Hamiltonian ***'
  ALLOCATE(mat_var(3)%non_zero_columns(n3d),                                            &
           mat_var(3)%row_index(lenbuf),                                                &
           mat_var(3)%mat_z%packed_columns(lenbuf),                                     &
           mat_var(3)%mat_z%matrix_diagonal(n3d) )
  drop_tol=drop_hamiltonian
  Call Pack_Triangle ('hamiltonian',                                                    &
                       mat_z%triangle_hamiltonian,                                      &
                       mat_var(3)%mat_z%packed_columns,                                 &
                       mat_var(3)%non_zero_columns,                                     &
                       mat_var(3)%row_index,                                            &
                       mat_var(3)%number,                                               &
                       mat_var(3)%mat_z%matrix_diagonal)
  DEALLOCATE(mat_var(3)%non_zero_columns,                                               &
             mat_var(3)%row_index,                                                      &
             mat_var(3)%mat_z%packed_columns,                                           &
             mat_var(3)%mat_z%matrix_diagonal )
!***********************************************************************
  END SUBROUTINE Read_Input_Matrices_z
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
  SUBROUTINE Diagonalize_Packed_Matrices_d (mat_d,matrix_source,mat_size,                    &
                                            get_eigenvectors,to_disk,                        &
                                            file_key,ham,s)
  IMPLICIT NONE
  TYPE(REAL_MATRICES)                       :: mat_d
  CHARACTER(LEN=*)                          :: matrix_source
  LOGICAL                                   :: get_eigenvectors 
  INTEGER                                   :: mat_size
  LOGICAL, OPTIONAL                         :: to_disk 
  CHARACTER (LEN=*), OPTIONAL               :: file_key 
  REAL(idp), DIMENSION(:), OPTIONAL         :: ham
  REAL(idp), DIMENSION(:), OPTIONAL         :: s
  REAL(idp), DIMENSION(:), ALLOCATABLE      :: eig
  REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: eigenvectors
  REAL(idp), DIMENSION(:), ALLOCATABLE      :: work
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
  ALLOCATE( eig(1:mat_size), work(1:lwork) )
  IF ( matrix_source /= 'input') THEN
       ALLOCATE( mat_d%triangle_hamiltonian(1:tri_size) )
       IF (non_orth) THEN
           ALLOCATE( mat_d%triangle_overlap(1:tri_size) )
       END IF
  END IF
  n_v='n'
  IF (get_eigenvectors == .true. ) THEN
      ALLOCATE(eigenvectors(1:mat_size,1:mat_size))
      n_v='v'
  END IF
  IF (matrix_source == 'disk' ) THEN
      write(iout,*)
      write(iout,*) '                    *** Diagonalizing Big Matrix From Disk ***'
      write(iout,*)
      IF(non_orth) THEN
         Call Fill_Matrix_from_Disk(mat_d%triangle_overlap,                                 &
                                    mat%non_zero_overlap_elements,                          &
                                    50)
         Call Fill_Matrix_from_Disk(mat_d%triangle_hamiltonian,                             &
                                    mat%non_zero_hamiltonian_elements,                      &
                                    50)
         call dspgv(1,n_v,'u',mat_size,mat_d%triangle_hamiltonian,mat_d%triangle_overlap,   &
                                  eig,eigenvectors,mat_size,work,info)
      ELSE 
         Call Fill_Matrix_from_Disk(mat_d%triangle_hamiltonian,                             &
                                    mat%non_zero_hamiltonian_elements,                      &
                                    50)
         call dspev(n_v,'u',mat_size,mat_d%triangle_hamiltonian,eig,eigenvectors,mat_size,work,info)
      END IF
  ELSE IF(matrix_source == 'buffers' ) THEN
      write(iout,*)
      write(iout,*) '                    *** Diagonalizing Big Matrix From Buffers ***'
      IF (non_orth) THEN
          Call IOsys('read integer number_non_zero_overlap_elements from '//                &
                     'packed_matrices',1,mat_var(1)%number,0,' ')
          ALLOCATE(mat_var(1)%non_zero_columns(mat_size),                                   &
                   mat_var(1)%row_index(mat_var(1)%number),                                 &
                   mat_var(1)%mat_d%packed_columns(mat_var(1)%number),                      &
                   mat_var(1)%mat_d%matrix_diagonal(mat_size) )
          Call Read_and_Write_Column_Packed_Matrices (                                      &
                                          mat_var(1)%mat_d%packed_columns,                  &
                                          mat_var(1)%non_zero_columns,                      &
                                          mat_var(1)%row_index,                             &
                                          'read',                                           &
                                          overlap,                                          &
                                          mat_var(1)%number,                                &
                                          mat_var(1)%mat_d%matrix_diagonal )
          Call Fill_Upper_Triangular_Matrix_from_Column_Buffers (                           &
                                          mat_d%triangle_overlap,                           &
                                          mat_var(1)%mat_d%matrix_diagonal,                 &
                                          mat_var(1)%mat_d%packed_columns,                  &
                                          mat_var(1)%non_zero_columns,                      &
                                          mat_var(1)%row_index)
          DEALLOCATE(mat_var(1)%row_index, mat_var(1)%mat_d%packed_columns )
          Call IOsys('read integer number_non_zero_hamiltonian_elements from '//            &
                     'packed_matrices',1,mat_var(1)%number,0,' ') 
          ALLOCATE(mat_var(1)%row_index(mat_var(1)%number),                                 &
                   mat_var(1)%mat_d%packed_columns(mat_var(1)%number) )
          Call Read_and_Write_Column_Packed_Matrices (                                      &
                                            mat_var(1)%mat_d%packed_columns,                &
                                            mat_var(1)%non_zero_columns,                    &
                                            mat_var(1)%row_index,                           &
                                            'read',                                         &
                                            hamiltonian,                                    &
                                            mat_var(1)%number,                              &
                                            mat_var(1)%mat_d%matrix_diagonal )
          Call Fill_Upper_Triangular_Matrix_from_Column_Buffers (                           &
                                              mat_d%triangle_hamiltonian,                   &
                                              mat_var(1)%mat_d%matrix_diagonal,             &
                                              mat_var(1)%mat_d%packed_columns,              &
                                              mat_var(1)%non_zero_columns,                  &
                                              mat_var(1)%row_index)
          call dspgv(1,n_v,'u',mat_size,mat_d%triangle_hamiltonian,                         &
                                  mat_d%triangle_overlap,                                   &
                                  eig,eigenvectors,mat_size,work,info)
      ELSE
          Call IOsys('read integer number_non_zero_hamiltonian_elements from '//            &
                     'packed_matrices',1,mat_var(1)%number,0,' ') 
          ALLOCATE(mat_var(1)%non_zero_columns(mat_size),                                   &
                   mat_var(1)%row_index(mat_var(1)%number),                                 &
                   mat_var(1)%mat_d%packed_columns(mat_var(1)%number),                      &
                   mat_var(1)%mat_d%matrix_diagonal(mat_size) )
          Call Read_and_Write_Column_Packed_Matrices (                                      &
                                            mat_var(1)%mat_d%packed_columns,                &
                                            mat_var(1)%non_zero_columns,                    &
                                            mat_var(1)%row_index,                           &
                                            'read',                                         &
                                            hamiltonian,                                    &
                                            mat_var(1)%number,                              &
                                            mat_var(1)%mat_d%matrix_diagonal )
          Call Fill_Upper_Triangular_Matrix_from_Column_Buffers (                           &
                                              mat_d%triangle_hamiltonian,                   &
                                              mat_var(1)%mat_d%matrix_diagonal,             &
                                              mat_var(1)%mat_d%packed_columns,              &
                                              mat_var(1)%non_zero_columns,                  &
                                              mat_var(1)%row_index)
          call dspev(n_v,'u',mat_size,mat_d%triangle_hamiltonian,eig,eigenvectors,          &
                     mat_size,work,info)
      END IF
      DEALLOCATE(mat_var(1)%non_zero_columns,                                               &
                 mat_var(1)%row_index,                                                      &
                 mat_var(1)%mat_d%packed_columns,                                           &
                 mat_var(1)%mat_d%matrix_diagonal )
  ELSE IF(matrix_source == 'input') THEN
      write(iout,*)
      write(iout,*) '                    *** Diagonalizing Big Matrix From Input Matrices ***'
      write(iout,*)
      IF (non_orth) THEN
          call dspgv(1,n_v,'u',mat_size,ham,s,eig,eigenvectors,mat_size,work,info)
      ELSE
          call dspev(n_v,'u',mat_size,ham,eig,eigenvectors,mat_size,work,info)
      END IF 
  END IF
  Call IOsys('write real "hamiltonian_eigenvalues_'//file_key//'" to '//'packed_matrices',  &
              mat_size,eig,0,' ')
  IF (get_eigenvectors == .true. ) THEN
      Call IOsys('write real "hamiltonian_eigenvectors_'//file_key//'" to '//'packed_matrices',  &
                  mat_size*mat_size,eigenvectors,0,' ')
      DEALLOCATE(eigenvectors)
  END IF
  Call Print_Matrix(type_real_vector,eig,title='eigenvalues')
  write(iout,*)
  write(iout,*) '          *** Done Diagonalizing Big Matrix for this Symmetry ***'
  write(iout,*)
  DEALLOCATE( eig, work )
  IF ( matrix_source /= 'input') THEN
       DEALLOCATE( mat_d%triangle_hamiltonian )
       IF (non_orth) THEN
           DEALLOCATE( mat_d%triangle_overlap )
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
  SUBROUTINE Diagonalize_Packed_Matrices_z (mat_z,matrix_source,mat_size,                    &
                                            get_eigenvectors,to_disk,                        &
                                            file_key,ham,s)
  IMPLICIT NONE
  TYPE(COMPLEX_MATRICES)                    :: mat_z
  CHARACTER(LEN=*)                          :: matrix_source
  LOGICAL                                   :: get_eigenvectors 
  INTEGER                                   :: mat_size
  LOGICAL, OPTIONAL                         :: to_disk 
  CHARACTER (LEN=*), OPTIONAL               :: file_key 
  COMPLEX(idp), DIMENSION(:), OPTIONAL      :: ham
  COMPLEX(idp), DIMENSION(:), OPTIONAL      :: s
  REAL(idp), DIMENSION(:), ALLOCATABLE      :: eig
  COMPLEX(idp), DIMENSION(:,:), ALLOCATABLE :: eigenvectors
  COMPLEX(idp), DIMENSION(:), ALLOCATABLE   :: work
  REAL(idp), DIMENSION(:), ALLOCATABLE      :: rwork
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
  ALLOCATE( eig(1:mat_size), work(1:lwork), rwork(1:lwork) )
  IF ( matrix_source /= 'input') THEN
       ALLOCATE( mat_z%triangle_hamiltonian(1:tri_size) )
       IF (non_orth) THEN
           ALLOCATE( mat_z%triangle_overlap(1:tri_size) )
       END IF
  END IF
  n_v='n'
  IF (get_eigenvectors == .true. ) THEN
      ALLOCATE(eigenvectors(1:mat_size,1:mat_size))
      n_v='v'
  END IF
  IF (matrix_source == 'disk' ) THEN
      write(iout,*)
      write(iout,*) '                    *** Diagonalizing Big Matrix From Disk ***'
      write(iout,*)
      IF(non_orth) THEN
         Call Fill_Matrix_from_Disk(mat_z%triangle_overlap,                                 &
                                    mat%non_zero_overlap_elements,                          &
                                    50)
         Call Fill_Matrix_from_Disk(mat_z%triangle_hamiltonian,                             &
                                    mat%non_zero_hamiltonian_elements,                      &
                                    50)
         call zhpgv(1,n_v,'u',mat_size,mat_z%triangle_hamiltonian,                          &
                                             mat_z%triangle_overlap,                        &
                                             eig,eigenvectors,mat_size,work,rwork,info)
      ELSE 
         Call Fill_Matrix_from_Disk(mat_z%triangle_hamiltonian,                             &
                                    mat%non_zero_hamiltonian_elements,                      &
                                    50)
         call zhpev(n_v,'u',mat_size,mat_z%triangle_hamiltonian,eig,eigenvectors,           &
                            mat_size,work,rwork,info)
      END IF
  ELSE IF(matrix_source == 'buffers' ) THEN
      write(iout,*)
      write(iout,*) '                    *** Diagonalizing Big Matrix From Buffers ***'
      IF (non_orth) THEN
          Call IOsys('read integer number_non_zero_overlap_elements from '//                &
                     'packed_matrices',1,mat_var(1)%number,0,' ')
          ALLOCATE(mat_var(1)%non_zero_columns(mat_size),                                   &
                   mat_var(1)%row_index(mat_var(1)%number),                                 &
                   mat_var(1)%mat_z%packed_columns(mat_var(1)%number),                      &
                   mat_var(1)%mat_z%matrix_diagonal(mat_size) )
          Call Read_and_Write_Column_Packed_Matrices (                                      &
                                          mat_var(1)%mat_z%packed_columns,                  &
                                          mat_var(1)%non_zero_columns,                      &
                                          mat_var(1)%row_index,                             &
                                          'read',                                           &
                                          overlap,                                          &
                                          mat_var(1)%number,                                &
                                          mat_var(1)%mat_z%matrix_diagonal )
          Call Fill_Upper_Triangular_Matrix_from_Column_Buffers (                           &
                                          mat_z%triangle_overlap,                           &
                                          mat_var(1)%mat_z%matrix_diagonal,                 &
                                          mat_var(1)%mat_z%packed_columns,                  &
                                          mat_var(1)%non_zero_columns,                      &
                                          mat_var(1)%row_index)
          DEALLOCATE(mat_var(1)%row_index, mat_var(1)%mat_z%packed_columns )
          Call IOsys('read integer number_non_zero_hamiltonian_elements from '//            &
                     'packed_matrices',1,mat_var(1)%number,0,' ') 
          ALLOCATE(mat_var(1)%row_index(mat_var(1)%number),                                 &
                   mat_var(1)%mat_z%packed_columns(mat_var(1)%number) )
          Call Read_and_Write_Column_Packed_Matrices (                                      &
                                            mat_var(1)%mat_z%packed_columns,                &
                                            mat_var(1)%non_zero_columns,                    &
                                            mat_var(1)%row_index,                           &
                                            'read',                                         &
                                            hamiltonian,                                    &
                                            mat_var(1)%number,                              &
                                            mat_var(1)%mat_z%matrix_diagonal )
          Call Fill_Upper_Triangular_Matrix_from_Column_Buffers (                           &
                                              mat_z%triangle_hamiltonian,                   &
                                              mat_var(1)%mat_z%matrix_diagonal,             &
                                              mat_var(1)%mat_z%packed_columns,              &
                                              mat_var(1)%non_zero_columns,                  &
                                              mat_var(1)%row_index) 
         call zhpgv(1,n_v,'u',mat_size,mat_z%triangle_hamiltonian,                          &
                                           mat_z%triangle_overlap,                          &
                                           eig,eigenvectors,mat_size,work,                  &
                                           rwork,info)
      ELSE
          Call IOsys('read integer number_non_zero_hamiltonian_elements from '//            &
                     'packed_matrices',1,mat_var(1)%number,0,' ') 
          ALLOCATE(mat_var(1)%non_zero_columns(mat_size),                                   &
                   mat_var(1)%row_index(mat_var(1)%number),                                 &
                   mat_var(1)%mat_z%packed_columns(mat_var(1)%number),                      &
                   mat_var(1)%mat_z%matrix_diagonal(mat_size) )
          Call Read_and_Write_Column_Packed_Matrices (                                      &
                                            mat_var(1)%mat_z%packed_columns,                &
                                            mat_var(1)%non_zero_columns,                    &
                                            mat_var(1)%row_index,                           &
                                            'read',                                         &
                                            hamiltonian,                                    &
                                            mat_var(1)%number,                              &
                                            mat_var(1)%mat_z%matrix_diagonal )
          Call Fill_Upper_Triangular_Matrix_from_Column_Buffers (                           &
                                              mat_z%triangle_hamiltonian,                   &
                                              mat_var(1)%mat_z%matrix_diagonal,             &
                                              mat_var(1)%mat_z%packed_columns,              &
                                              mat_var(1)%non_zero_columns,                  &
                                              mat_var(1)%row_index)
          Call zhpev(n_v,'u',mat_size,mat_z%triangle_hamiltonian,                           &
                                      eig,eigenvectors,mat_size,work,rwork,info)
      END IF
      DEALLOCATE(mat_var(1)%non_zero_columns,                                               &
                 mat_var(1)%row_index,                                                      &
                 mat_var(1)%mat_z%packed_columns,                                           &
                 mat_var(1)%mat_z%matrix_diagonal )
  ELSE IF(matrix_source == 'input') THEN
      write(iout,*)
      write(iout,*) '                    *** Diagonalizing Big Matrix From Input Matrices ***'
      write(iout,*)
      IF (non_orth) THEN
          call zhpgv(1,n_v,'u',mat_size,ham,s,eig,eigenvectors,mat_size,work,rwork,info)
      ELSE
          call zhpev(n_v,'u',mat_size,ham,eig,eigenvectors,mat_size,work,rwork,info)
      END IF 
  END IF
  Call IOsys('write real "hamiltonian_eigenvalues_'//file_key//'" to '//'packed_matrices',  &
              mat_size,eig,0,' ')
  IF (get_eigenvectors == .true. ) THEN
      Call IOsys('write real "hamiltonian_eigenvectors_'//file_key//'" to '//'packed_matrices',  &
                  2*mat_size*mat_size,eigenvectors,0,' ')
      DEALLOCATE(eigenvectors)
  END IF
  Call Print_Matrix(type_real_vector,eig,title='eigenvalues')
  write(iout,*)
  write(iout,*) '          *** Done Diagonalizing Big Matrix for this Symmetry ***'
  write(iout,*)
  DEALLOCATE( eig, work, rwork )
  IF ( matrix_source /= 'input') THEN
       DEALLOCATE( mat_z%triangle_hamiltonian )
       IF (non_orth) THEN
           DEALLOCATE( mat_z%triangle_overlap )
       END IF
  END IF
!***********************************************************************
  END SUBROUTINE Diagonalize_Packed_Matrices_z
!***********************************************************************
!***********************************************************************
  END  MODULE Packed_Matrix_Module
!***********************************************************************
!***********************************************************************
