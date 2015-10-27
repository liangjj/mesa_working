
!***********************************************************************
                           MODULE Matrix_Utility_Module
                            USE input_output
                            USE dvrprop_global,                                          &
                                ONLY : n3d
                            USE Global_Time_Propagation_Module,                          &
                                ONLY : print_cc, print_packed_matrices, print_buffers,   &
                                       print_input_matrices, print_internal_matrices,    &                                      
                                       print_packed_disk_buffers 
                            IMPLICIT NONE
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
                          INTERFACE Read_and_Write_Column_Packed_Matrices
                   MODULE PROCEDURE Read_and_Write_Column_Packed_Matrices_d,             &
                                    Read_and_Write_Column_Packed_Matrices_z         
                          END INTERFACE Read_and_Write_Column_Packed_Matrices
!
                          INTERFACE Read_and_Write_Packed_Matrices
                   MODULE PROCEDURE Read_and_Write_Packed_Matrices_d,                    &
                                    Read_and_Write_Packed_Matrices_z
                          END INTERFACE Read_and_Write_Packed_Matrices
!
                          INTERFACE Fill_Matrix_from_Buffers
                   MODULE PROCEDURE Fill_Triangular_Matrix_from_Buffers_d,               &
                                    Fill_Triangular_Matrix_from_Buffers_z,               &              
                                    Fill_Square_Matrix_from_Buffers_d,                   &
                                    Fill_Square_Matrix_from_Buffers_z                       
                          END INTERFACE Fill_Matrix_from_Buffers
!
                          INTERFACE Fill_Upper_Triangular_Matrix_from_Column_Buffers
                   MODULE PROCEDURE Fill_Upper_Triangular_Matrix_from_Column_Buffers_d,  &
                                    Fill_Upper_Triangular_Matrix_from_Column_Buffers_z                 
                          END INTERFACE Fill_Upper_Triangular_Matrix_from_Column_Buffers
!
                          INTERFACE Fill_Lower_Triangular_Matrix_from_Column_Buffers
                   MODULE PROCEDURE Fill_Lower_Triangular_Matrix_from_Column_Buffers_d,  &
                                    Fill_Upper_Triangular_Matrix_from_Column_Buffers_z                 
                          END INTERFACE Fill_Lower_Triangular_Matrix_from_Column_Buffers
!
                          INTERFACE Fill_General_Matrix_from_Column_Buffers
                   MODULE PROCEDURE Fill_General_Matrix_from_Column_Buffers_d,           &
                                    Fill_General_Matrix_from_Column_Buffers_z                 
                          END INTERFACE Fill_General_Matrix_from_Column_Buffers
!
                          INTERFACE Fill_Buffers_from_Matrix
                   MODULE PROCEDURE Fill_Buffers_from_Triangular_Matrix_d,               &
                                    Fill_Buffers_from_Triangular_Matrix_z
                          END INTERFACE Fill_Buffers_from_Matrix
!
                          INTERFACE Fill_Buffers_from_Disk
                   MODULE PROCEDURE Fill_Buffers_from_Disk_d,                            &
                                    Fill_Buffers_from_Disk_z
                          END INTERFACE Fill_Buffers_from_Disk
!
                          INTERFACE Fill_Matrix_from_Disk
                   MODULE PROCEDURE Fill_Square_Matrix_from_Disk_d,                      &
                                    Fill_Square_Matrix_from_Disk_z,                      &
                                    Fill_Triangular_Matrix_from_Disk_d,                  &    
                                    Fill_Triangular_Matrix_from_Disk_z                 
                          END INTERFACE Fill_Matrix_from_Disk
!
                          INTERFACE Fill_Matrix_from_Input
                   MODULE PROCEDURE Fill_Square_Matrix_from_Input_d,                     &
                                    Fill_Square_Matrix_from_Input_z,                     &
                                    Fill_Triangular_Matrix_from_Input_d,                 &    
                                    Fill_Triangular_Matrix_from_Input_z                 
                          END INTERFACE Fill_Matrix_from_Input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           Contains
!***********************************************************************
!***********************************************************************
!***begin prologue     Read_and_Write_Column_Packed_Matrices_d
!***date written       000710   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           one-dim
!***author             schneider, barry (nsf)
!***source             
!***purpose            Create or open and read column packed files. 
!***                   What is stored are the non zero elements of each column
!***                   of the upper triangle.  The array packed_columns has
!***                   these.  The row index is stored is row_index and the
!***                   number of non zero elements in the column in non_zero_column.
!***                   Storage order of elements is a11, a21, a31, ....an1 etc.
!***                   if the entire column is stored.  In some cases only the
!***                   triangle is stored and often without diagonal.  Then its 
!***                   a11....an1,a22,....an2,a33... or with the diagonal missing.
!***references         
!
!***routines called    
!***end prologue       Read_and_Write_Column_Packed_Matrices_d
  Subroutine Read_and_Write_Column_Packed_Matrices_d (                            &
                                           packed_columns,                        &
                                           non_zero_columns,                      &
                                           row_index,                             &
                                           read_or_write,                         &
                                           m_type,                                &
                                           number,                                &
                                           matrix_diagonal)
  IMPLICIT NONE
  REAL*8,   DIMENSION(:), OPTIONAL  :: matrix_diagonal
  REAL*8,   DIMENSION(:)            :: packed_columns
  INTEGER,  DIMENSION(:)            :: non_zero_columns
  INTEGER,  DIMENSION(:)            :: row_index
  INTEGER                           :: number
  INTEGER                           :: column_length
  CHARACTER(LEN=*)                  :: m_type
  CHARACTER(LEN=*)                  :: read_or_write
!
  column_length=size(non_zero_columns)
  IF (read_or_write == 'write') Then
      Call IOsys('write integer "number_non_zero_'//m_type//'_elements" to '//   &
                 'packed_matrices',1,number,0,' ')
      write(iout,*) ' Writing packed file = '//m_type
      write(iout,*) ' Number of  Non Zero Elements = ',number
      Call IOsys('write integer "non_zero_'//m_type//'_columns" to '//           &
                 'packed_matrices',column_length,non_zero_columns,0,' ')
      Call IOsys('write integer "'//m_type//'_row_index" to '//                  &
                 'packed_matrices',number,row_index,0,' ')
      Call IOsys('write real "'//m_type//'_packed_columns" to '//                &
                 'packed_matrices',number,packed_columns,0,' ')
      IF ( PRESENT (matrix_diagonal) ) THEN
           Call IOsys('write real "'//m_type//'_matrix_diagonal" to '//          &
                      'packed_matrices',column_length,matrix_diagonal,0,' ')
      END IF
  ELSE IF (read_or_write == 'read') THEN
!      Call IOsys('read integer "number_non_zero_'//m_type//'_elements" from '//  &
!                 'packed_matrices',1,number,0,' ')
      write(iout,*) ' Reading packed file = '//m_type
      write(iout,*) ' Number of  Non Zero Elements = ',number
      Call IOsys('read integer "non_zero_'//m_type//'_columns" from '//           &
                'packed_matrices',column_length,non_zero_columns,0,' ')
      Call IOsys('read integer "'//m_type//'_row_index" from '//                 &
                 'packed_matrices',number,row_index,0,' ')
      Call IOsys('read real "'//m_type//'_packed_columns" from '//               &
                 'packed_matrices',number,packed_columns,0,' ')
      IF ( PRESENT (matrix_diagonal) ) THEN
           Call IOsys('read real "'//m_type//'_matrix_diagonal" from '//         &
                      'packed_matrices',column_length,matrix_diagonal,0,' ')
      END IF
  END IF
  IF(print_packed_disk_buffers) THEN
     write(iout,*) 'number of non zero elements = ', number
     write(iout,*) 'non zero columns = ', non_zero_columns
     write(iout,*) 'row index = ', row_index
     write(iout,*) 'packed columns = ', packed_columns
     write(iout,*) 'diagonals = ', matrix_diagonal
  END IF
END SUBROUTINE Read_and_Write_Column_Packed_Matrices_d 
!***********************************************************************
!***********************************************************************
!***begin prologue     Read_and_Write_Column_Packed_Matrices_z 
!***date written       000710   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           one-dim
!***author             schneider, barry (nsf)
!***source             
!***purpose            Create or open and read column packed files. 
!***                   What is stored are the non zero elements of each column
!***                   of the upper triangle.  The array packed_columns has
!***                   these.  The row index is stored is row_index and the
!***                   number of non zero elements in the column in non_zero_column.
!***                   Storage is a11, a12, a22, a13, a23, a33 etc.
!***references         
!
!***routines called    
!***end prologue      Read_and_Write_Column_Packed_Matrices_z  
  Subroutine Read_and_Write_Column_Packed_Matrices_z (                         &
                                           packed_columns,                     &
                                           non_zero_columns,                   &
                                           row_index,                          &
                                           read_or_write,                      &
                                           m_type,                             &
                                           number,                             &
                                           matrix_diagonal)
  IMPLICIT NONE
  COMPLEX*16,   DIMENSION(:), OPTIONAL   :: matrix_diagonal
  COMPLEX*16,   DIMENSION(:)             :: packed_columns
  INTEGER,      DIMENSION(:)             :: non_zero_columns
  INTEGER,      DIMENSION(:)             :: row_index
  INTEGER                                :: number
  INTEGER                                :: column_length
  CHARACTER(LEN=*)                       :: m_type
  CHARACTER(LEN=*)                       :: read_or_write
!
  column_length=size(non_zero_columns)
  IF (read_or_write == 'write') Then
      Call IOsys('write integer "number_non_zero_'//m_type//'_elements" to '//   &
                 'packed_matrices',1,number,0,' ')
      write(iout,*) ' Writing packed file = '//m_type
      write(iout,*) ' Number of  Non Zero Elements = ',number
      Call IOsys('write integer "non_zero_'//m_type//'_columns" to '//           &
                 'packed_matrices',column_length,non_zero_columns,0,' ')
      Call IOsys('write integer "'//m_type//'_row_index" to '//                  &
                 'packed_matrices',number,row_index,0,' ')
      Call IOsys('write real "'//m_type//'_packed_columns" to '//                &
                 'packed_matrices',2*number,packed_columns,0,' ')
      IF ( PRESENT (matrix_diagonal) ) THEN
           Call IOsys('write real "'//m_type//'_matrix_diagonal" to '//          & 
                      'packed_matrices',2*column_length,matrix_diagonal,0,' ')
      END IF
  ELSE IF (read_or_write == 'read') THEN
!      Call IOsys('read integer "number_non_zero_'//m_type//'_elements" from '//  &
!                 'packed_matrices',1,number,0,' ')
      write(iout,*) ' Reading packed file = '//m_type
      write(iout,*) ' Number of  Non Zero Elements = ',number
      Call IOsys('read integer "non_zero_'//m_type//'_columns" from '//          &
                 'packed_matrices',column_length,non_zero_columns,0,' ')
      Call IOsys('read integer "'//m_type//'_row_index" from '//                 &
                 'packed_matrices',number,row_index,0,' ')
      Call IOsys('read real "'//m_type//'_packed_columns" from '//               &
                 'packed_matrices',2*number,packed_columns,0,' ')
      IF ( PRESENT (matrix_diagonal) ) THEN
           Call IOsys('read real "'//m_type//'_matrix_diagonal" from '//         &
                      'packed_matrices',2*column_length,matrix_diagonal,0,' ')
      END IF
  END IF
  IF(print_packed_disk_buffers) THEN
     write(iout,*) 'number of non zero elements = ', number
     write(iout,*) 'non zero columns = ', non_zero_columns
     write(iout,*) 'row index = ', row_index
     write(iout,*) 'packed columns = ', packed_columns
     write(iout,*) 'diagonals = ', matrix_diagonal
  END IF
END SUBROUTINE Read_and_Write_Column_Packed_Matrices_z 
!***********************************************************************
!***********************************************************************
!Deck Read_and_Write_Packed_Matrices_d
!***begin prologue     Read_and_Write_Packed_Matrices_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Read_and_Write packed files stored in the old two index
!***                   method with separate storage of diagonals.
!***references
!***routines called
!***end prologue       Read_and_Write_Packed_Matrices_d
!
  SUBROUTINE Read_and_Write_Packed_Matrices_d(ibuf,                             &
                                              h_buf,                            &
                                              read_or_write,                    &
                                              m_type,                           &
                                              non_zero,                         &
                                              matrix_diagonal )
  IMPLICIT NONE
  INTEGER, DIMENSION(2,*)               :: ibuf
  REAL*8, DIMENSION(:)                  :: h_buf
  REAL*8, DIMENSION(:), OPTIONAL        :: matrix_diagonal
  CHARACTER(LEN=*)                      :: m_type 
  CHARACTER(LEN=*)                      :: read_or_write 
  INTEGER                               :: non_zero
  IF (read_or_write == 'read') THEN
!      Call IOsys('read integer "number_non_zero_'//m_type//'_elements" from '//   &
!                 'packed_matrices',1,non_zero,0,' ')
      Call IOsys('read integer "'//m_type//'_integer_buffers" from '//            &
                 'packed_matrices',2*non_zero,ibuf,0,' ') 
      Call IOsys('read real "'//m_type//'_real_buffers" from '//                  &
                 'packed_matrices',non_zero,h_buf,0,' ') 
      IF ( PRESENT (matrix_diagonal) ) THEN
           Call IOsys('read real "'//m_type//'_diagonals" from '//                &
                      'packed_matrices',n3d,matrix_diagonal,0,' ') 
      END IF
  ELSE
      Call IOsys('write integer "number_non_zero_'//m_type//'_elements" to '//    &
                 'packed_matrices',1,non_zero,0,' ') 
      Call IOsys('write integer "'//m_type//'_integer_buffers" to '//             &
                 'packed_matrices',2*non_zero,ibuf,0,' ') 
      Call IOsys('write real "'//m_type//'_real_buffers" to '//                   &
                 'packed_matrices',non_zero,h_buf,0,' ') 
      IF ( PRESENT (matrix_diagonal) ) THEN
           Call IOsys('write real "'//m_type//'_diagonals" to '//                 &
                      'packed_matrices',n3d,matrix_diagonal,0,' ') 
      END IF
  END IF
!***********************************************************************
  END SUBROUTINE Read_and_Write_Packed_Matrices_d
!***********************************************************************
!***********************************************************************
!Deck Read_and_Write_Packed_Matrices_z
!***begin prologue     Read_and_Write_Packed_Matrices_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       Read_and_Write_Packed_Matrices_z
!
  SUBROUTINE Read_and_Write_Packed_Matrices_z(ibuf,                             &
                                              h_buf,                            &
                                              read_or_write,                    &
                                              m_type,                           &
                                              non_zero,                         &
                                              matrix_diagonal )
  IMPLICIT NONE
  INTEGER, DIMENSION(2,*)               :: ibuf
  COMPLEX*16, DIMENSION(:)              :: h_buf
  COMPLEX*16, DIMENSION(:), OPTIONAL    :: matrix_diagonal
  CHARACTER(LEN=*)                      :: m_type 
  CHARACTER(LEN=*)                      :: read_or_write 
  INTEGER                               :: non_zero
  IF (read_or_write == 'read') THEN
!      Call IOsys('read integer "non_zero_'//m_type//'_elements" from '//          &
!                 'packed_matrices',1,non_zero,0,' ')
      Call IOsys('read integer "'//m_type//'_integer_buffers" from '//            &
                 'packed_matrices',2*non_zero,ibuf,0,' ') 
      Call IOsys('read real "'//m_type//'_real_buffers" from '//                  &
                 'packed_matrices',2*non_zero,h_buf,0,' ') 
      IF ( PRESENT (matrix_diagonal) ) THEN
           Call IOsys('read real "'//m_type//'_diagonals" from '//                &
                      'packed_matrices',2*n3d,matrix_diagonal,0,' ') 
      END IF
  ELSE
      Call IOsys('write integer "non_zero_'//m_type//'_elements" to '//           &
                 'packed_matrices',1,non_zero,0,' ') 
      Call IOsys('write integer "'//m_type//'_integer_buffers" to '//             &
                 'packed_matrices',2*non_zero,ibuf,0,' ') 
      Call IOsys('write real "'//m_type//'_real_buffers" to '//                   &
                 'packed_matrices',2*non_zero,h_buf,0,' ') 
      IF ( PRESENT (matrix_diagonal) ) THEN
           Call IOsys('write real "'//m_type//'_diagonals" to '//                 &
                 'packed_matrices',2*n3d,matrix_diagonal,0,' ') 
      END IF 
  END IF
! ***********************************************************************
  END SUBROUTINE Read_and_Write_Packed_Matrices_z
!***********************************************************************
!***********************************************************************
!Deck Fill_Triangular_Matrix_from_Buffers_d
!***begin prologue     Fill_Triangular_Matrix_from_Buffers_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill triangular matrix from a buffered array.
!***                   The buffers contain the packed matrix in the
!***                   two non-zero index format.
!***references
!***routines called
!***end prologue       Fill_Triangular_Matrix_d
!
  SUBROUTINE Fill_Triangular_Matrix_from_Buffers_d(matrix,ibuf,h_buf,          &
                                                   matrix_diagonal,non_zero)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                  :: matrix
  INTEGER, DIMENSION(2,*)               :: ibuf
  REAL*8, DIMENSION(:)                  :: h_buf
  REAL*8, DIMENSION(:)                  :: matrix_diagonal
  INTEGER                               :: non_zero
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: ind_i
  INTEGER                               :: ind_j
  INTEGER                               :: ij
  CHARACTER(LEN=80)                     :: title
  matrix(:) = 0.d0
  DO i=1,non_zero
     ind_i = ibuf(1,i)
     ind_j = ibuf(2,i)
     ij = ind_i * (ind_i - 1 ) / 2 + ind_j
     matrix(ij) = h_buf(i)
  END DO
  ij = 0
  DO i = 1, n3d
     ij = ij + i
     matrix(ij) = matrix_diagonal(i)
  END DO 
  IF (print_internal_matrices) THEN
      title= 'Printing Triangular Matrix from Buffers Subroutine'     
      Write(iout,*) title
      ij = 0
      DO i= 1, n3d
         write(iout,1) i
         write(iout,2) ( matrix(j), j=ij+1, ij + i )
         ij  = ij  + i
      END DO
  END IF
1 Format(/,5x,'Row = ',i4)
2 Format( (15x,5e15.8) )
!************* **********************************************************
  END SUBROUTINE Fill_Triangular_Matrix_from_Buffers_d
!***********************************************************************
!***********************************************************************
!Deck Fill_Triangular_Matrix_from_Buffers_z
!***begin prologue     Fill_Triangular_Matrix_from_Buffers_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill_Triangular
!***
!***references
!***routines called
!***end prologue       Fill_Triangular_Matrix_z
!
  SUBROUTINE Fill_Triangular_Matrix_from_Buffers_z(matrix,ibuf,h_buf,matrix_diagonal,     &
                                                   non_zero)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)              :: matrix
  INTEGER, DIMENSION(2,*)               :: ibuf
  COMPLEX*16, DIMENSION(:)              :: h_buf
  COMPLEX*16, DIMENSION(:)              :: matrix_diagonal
  INTEGER                               :: non_zero
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: ind_i
  INTEGER                               :: ind_j
  INTEGER                               :: ij
  CHARACTER(LEN=80)                     :: title
  matrix(:) = 0.d0
  DO i=1,non_zero
     ind_i = ibuf(1,i)
     ind_j = ibuf(2,i)
     ij = ind_i * (ind_i - 1 ) / 2 + ind_j
     matrix(ij) = h_buf(i)
  END DO
  ij = 0
  DO i = 1, n3d
     ij = ij + i
     matrix(ij) = matrix_diagonal(i)
  END DO 
  IF (print_internal_matrices) THEN
      title='Printing Matrix from Fill_Triangular_Matrix_from_Buffers Subroutine'     
      write(iout,*) title
      ij = 0
      DO i= 1, n3d
         write(iout,1) i
         write(iout,2) ( matrix(j), j=ij+1, ij + i )
         ij  = ij  + i
      END DO
  END IF
1 Format(/,5x,'Row = ',i4)
2 Format( (15x,4e15.8) )
!***********************************************************************
  END SUBROUTINE Fill_Triangular_Matrix_from_Buffers_z
!***********************************************************************
!***********************************************************************
!Deck Fill_Upper_Triangular_Matrix_from_Column_Buffers_d
!***begin prologue     Fill_Upper_Triangular_Matrix_from_Column_Buffers_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill upper triangular matrix from a column buffered
!***                   matrix.
!***references
!***routines called
!***end prologue       Fill_Upper_Triangular_Matrix_d
!
  Subroutine Fill_Upper_Triangular_Matrix_from_Column_Buffers_d (                 &
                                           upper_matrix,                          &
                                           matrix_diagonal,                       &
                                           packed_columns,                        &
                                           non_zero_columns,                      &
                                           row_index )
  IMPLICIT NONE
  REAL*8,   DIMENSION(:)                 :: upper_matrix
  REAL*8,   DIMENSION(:)                 :: matrix_diagonal
  REAL*8,   DIMENSION(:)                 :: packed_columns
  INTEGER,  DIMENSION(:)                 :: non_zero_columns
  INTEGER,  DIMENSION(:)                 :: row_index
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: jj
  INTEGER                                :: mat_ind
  INTEGER                                :: count
  CHARACTER(LEN=80)                      :: title
  upper_matrix(:) = 0.d0
  count = 0  
  DO j = 1 , n3d
     jj= j * ( j - 1 ) / 2
     DO i = 1 , non_zero_columns(j)
        count = count + 1
        mat_ind = row_index(count) + jj
        upper_matrix(mat_ind) = packed_columns(count)
     END DO
     upper_matrix(jj + j) = matrix_diagonal(j)
  END DO
  IF (print_internal_matrices) THEN
      title= 'Printing Triangular Matrix from Column_Buffers Subroutine'     
      Write(iout,*) title
      jj = 0
      DO i= 1, n3d
         write(iout,1) i
         write(iout,2) ( upper_matrix(j), j=jj + 1, jj + i )
         jj  = jj  + i
      END DO
  END IF
1 Format(/,5x,'Col = ',i4)
2 Format( (15x,5e15.8) )
!************* **********************************************************
  END SUBROUTINE Fill_Upper_Triangular_Matrix_from_Column_Buffers_d
!***********************************************************************
!***********************************************************************
!Deck Fill_Upper_Triangular_Matrix_from_Column_Buffers_z
!***begin prologue     Fill_Upper_Triangular_Matrix_from_Column_Buffers_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill_Triangular_Matrix
!***
!***references
!***routines called
!***end prologue       Fill_Upper_Triangular_Matrix_z
!
  Subroutine Fill_Upper_Triangular_Matrix_from_Column_Buffers_z (                 &
                                           upper_matrix,                          &
                                           matrix_diagonal,                       &
                                           packed_columns,                        &
                                           non_zero_columns,                      &
                                           row_index )
  IMPLICIT NONE
  COMPLEX*16,   DIMENSION(:)             :: upper_matrix
  COMPLEX*16,   DIMENSION(:)             :: matrix_diagonal
  COMPLEX*16,   DIMENSION(:)             :: packed_columns
  INTEGER,  DIMENSION(:)                 :: non_zero_columns
  INTEGER,  DIMENSION(:)                 :: row_index
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: jj
  INTEGER                                :: mat_ind
  INTEGER                                :: count
  CHARACTER(LEN=80)                      :: title
  upper_matrix(:) = 0.d0
  count = 0  
  DO j = 1 , n3d
     jj = j * ( j - 1 ) / 2
     DO i = 1 , non_zero_columns(j)
        count = count + 1
        mat_ind = row_index(count) + jj
        upper_matrix(mat_ind) = packed_columns(count)
     END DO
     upper_matrix(jj + j) = matrix_diagonal(j)
  END DO
  IF (print_internal_matrices) THEN
      title= 'Printing Triangular Matrix from Column_Buffers Subroutine'     
      Write(iout,*) title
      jj = 0
      DO i= 1, n3d
         write(iout,1) i
         write(iout,2) ( upper_matrix(j), j=jj + 1, jj + i )
         jj  = jj  + i
      END DO
  END IF
1 Format(/,5x,'Row = ',i4)
2 Format( (15x,6e15.8) )
!***********************************************************************
  END SUBROUTINE Fill_Upper_Triangular_Matrix_from_Column_Buffers_z
!***********************************************************************
!***********************************************************************
!Deck Fill_Lower_Triangular_Matrix_from_Column_Buffers_d
!***begin prologue     Fill_Lower_Triangular_Matrix_from_Column_Buffers_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill lower triangular matrix from a column packed matrix.
!***
!***references
!***routines called
!***end prologue       Fill_Lower_Triangular_Matrix_d
!
  Subroutine Fill_Lower_Triangular_Matrix_from_Column_Buffers_d (                 &
                                           lower_matrix,                          &
                                           matrix_diagonal,                       &
                                           packed_columns,                        &
                                           non_zero_columns,                      &
                                           row_index )
  IMPLICIT NONE
  REAL*8,   DIMENSION(:)                 :: lower_matrix
  REAL*8,   DIMENSION(:)                 :: matrix_diagonal
  REAL*8,   DIMENSION(:)                 :: packed_columns
  INTEGER,  DIMENSION(:)                 :: non_zero_columns
  INTEGER,  DIMENSION(:)                 :: row_index
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: jj
  INTEGER                                :: j_ind
  INTEGER                                :: mat_ind
  INTEGER                                :: count
  CHARACTER(LEN=80)                      :: title
  lower_matrix(:) = 0.d0
  count = 0  
  DO j = 1 , n3d
     jj = ( j - 1 ) * ( n3d + n3d - j) / 2
     DO i = 1 , non_zero_columns(j)
        count = count + 1
        mat_ind = row_index(count) + jj
        lower_matrix(mat_ind) = packed_columns(count)
     END DO
     lower_matrix(jj+j) = matrix_diagonal(j)
  END DO
  IF (print_internal_matrices) THEN
      title= 'Printing Triangular Matrix from Column_Buffers Subroutine'     
      Write(iout,*) title
      jj = 0
      DO i= 1, n3d
         write(iout,1) i
         write(iout,2) ( lower_matrix(j), j=jj + 1, jj + i )
         jj  = jj  + i
      END DO
  END IF
1 Format(/,5x,'Row = ',i4)
2 Format( (15x,5e15.8) )
!************* **********************************************************
  END SUBROUTINE Fill_Lower_Triangular_Matrix_from_Column_Buffers_d
!***********************************************************************
!***********************************************************************
!Deck Fill_Lower_Triangular_Matrix_from_Column_Buffers_z
!***begin prologue     Fill_Lower_Triangular_Matrix_from_Column_Buffers_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill_Triangular_Matrix
!***
!***references
!***routines called
!***end prologue       Fill_Lower_Triangular_Matrix_z
!
  Subroutine Fill_Lower_Triangular_Matrix_from_Column_Buffers_z (                 &
                                           matrix,                                &
                                           matrix_diagonal,                       &
                                           packed_columns,                        &
                                           non_zero_columns,                      &
                                           row_index )
  IMPLICIT NONE
  COMPLEX*16,   DIMENSION(:)             :: matrix
  COMPLEX*16,   DIMENSION(:)             :: matrix_diagonal
  COMPLEX*16,   DIMENSION(:)             :: packed_columns
  INTEGER,  DIMENSION(:)                 :: non_zero_columns
  INTEGER,  DIMENSION(:)                 :: row_index
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: jj
  INTEGER                                :: mat_ind
  INTEGER                                :: count
  CHARACTER(LEN=80)                      :: title
  matrix(:) = 0.d0
  count = 0  
  DO j = 1 , n3d
     jj = ( j - 1 ) * ( n3d + n3d - j) / 2
     DO i = 1 , non_zero_columns(j)
        count = count + 1
        mat_ind = row_index(count) + jj
        matrix(mat_ind) = packed_columns(count)
     END DO
     matrix(jj+j) = matrix_diagonal(j)
  END DO
  IF (print_internal_matrices) THEN
      title= 'Printing Triangular Matrix from Column_Buffers Subroutine'     
      Write(iout,*) title
      jj = 0
      DO i= 1, n3d
         write(iout,1) i
         write(iout,2) ( matrix(j), j=jj + 1, jj + i )
         jj  = jj  + i
      END DO
  END IF
1 Format(/,5x,'Row = ',i4)
2 Format( (15x,6e15.8) )
!***********************************************************************
  END SUBROUTINE Fill_Lower_Triangular_Matrix_from_Column_Buffers_z
!***********************************************************************
!***********************************************************************
!Deck Fill_General_Matrix_from_Column_Buffers_d
!***begin prologue     Fill_General_Matrix_from_Column_Buffers_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill General Matrix from Column Buffers
!***
!***references
!***routines called
!***end prologue       Fill_General_Matrix_from_Column_Buffers_d
!
  Subroutine Fill_General_Matrix_from_Column_Buffers_d (                    &
                                                        matrix,             &
                                                        packed_columns,     &
                                                        non_zero_columns,   &
                                                        row_index )
  IMPLICIT NONE
  REAL*8,   DIMENSION(:,:)                 :: matrix
  REAL*8,   DIMENSION(:)                   :: packed_columns
  INTEGER,  DIMENSION(:)                   :: non_zero_columns
  INTEGER,  DIMENSION(:)                   :: row_index
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: n_col
  INTEGER                                  :: n_row
  INTEGER                                  :: count
  CHARACTER(LEN=80)                        :: title
  matrix(:,:) = 0.d0
  n_col=size(matrix,2)
  n_row=size(matrix,1)
  count = 0  
  DO j = 1 , n_col
     DO i = 1 , non_zero_columns(j)
        count = count + 1
        matrix(row_index(count),j) = packed_columns(count)
     END DO
  END DO
  IF (print_internal_matrices) THEN
      title= 'Printing Matrix from Column_Buffers'     
      call prntfmn(title,matrix,n_col,n_row,n_col,n_row,iout,'e')
  END IF
!************* **********************************************************
  END SUBROUTINE Fill_General_Matrix_from_Column_Buffers_d
!***********************************************************************
!***********************************************************************
!Deck Fill_General_Matrix_from_Column_Buffers_z
!***begin prologue     Fill_General_Matrix_from_Column_Buffers_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill General Matrix from Column Buffers
!***
!***references
!***routines called
!***end prologue       Fill_General_Matrix_from_Column_Buffers_z
!
  Subroutine Fill_General_Matrix_from_Column_Buffers_z (                    &
                                                        matrix,             &
                                                        packed_columns,     &
                                                        non_zero_columns,   &
                                                        row_index )
  IMPLICIT NONE
  COMPLEX*16,   DIMENSION(:,:)                 :: matrix
  COMPLEX*16,   DIMENSION(:)                   :: packed_columns
  INTEGER,      DIMENSION(:)                   :: non_zero_columns
  INTEGER,      DIMENSION(:)                   :: row_index
  INTEGER                                      :: i
  INTEGER                                      :: j
  INTEGER                                      :: n_col
  INTEGER                                      :: n_row
  INTEGER                                      :: count
  CHARACTER(LEN=80)                            :: title
  matrix(:,:) = 0.d0
  n_col=size(matrix,2)
  n_row=size(matrix,1)
  count = 0  
  DO j = 1 , n_col
     DO i = 1 , non_zero_columns(j)
        count = count + 1
        matrix(row_index(count),j) = packed_columns(count)
     END DO
  END DO
  IF (print_internal_matrices) THEN
      title= 'Printing Matrix from Column_Buffers'     
      call prntcmn(title,matrix,n_col,n_row,n_col,n_row,iout,'e')
  END IF
!************************************************************************
  END SUBROUTINE Fill_General_Matrix_from_Column_Buffers_z
!***********************************************************************
!***********************************************************************
!Deck Fill_Square_Matrix_from_Buffers_d
!***begin prologue     Fill_Square_Matrix_from_Buffers_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill_Square_Matrix_from_Buffers using the two 
!***                   non zero index method of storage.
!***references
!***routines called
!***end prologue       Fill_Square_Matrix_from_Buffers_d
!
  SUBROUTINE Fill_Square_Matrix_from_Buffers_d(matrix,ibuf,h_buf,matrix_diagonal,      &
                                               non_zero )
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)                :: matrix    
  INTEGER, DIMENSION(2,*)               :: ibuf
  REAL*8, DIMENSION(:)                  :: h_buf
  REAL*8, DIMENSION(:)                  :: matrix_diagonal
  INTEGER                               :: non_zero
  INTEGER                               :: i
  INTEGER                               :: ind_i
  INTEGER                               :: ind_j
  CHARACTER(LEN=80)                     :: title
  matrix(:,:) = 0.d0
  DO i=1,non_zero
     ind_i = ibuf(1,i)
     ind_j = ibuf(2,i)
     matrix(ind_i,ind_j) = h_buf(i)
     matrix(ind_j,ind_i) = h_buf(i)
  END DO
  DO i = 1,n3d
     matrix(i,i) = matrix_diagonal(i)
  END DO
  IF (print_internal_matrices) THEN
      title= 'Printing Matrix from Fill_Square_Matrix_from_Buffers Subroutine'     
      Call prntfmn(title,matrix,n3d,n3d,n3d,n3d,iout,'e')
  END IF
!***********************************************************************
  END SUBROUTINE Fill_Square_Matrix_from_Buffers_d
!***********************************************************************
!***********************************************************************
!Deck Fill_Square_Matrix_from_Buffers_z
!***begin prologue     Fill_Square_Matrix_from_Buffers_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill_Square_Matrix_from_Buffers_z
!***
!***references
!***routines called
!***end prologue       Fill_Square_Matrix_from_Buffers_z
!
  SUBROUTINE Fill_Square_Matrix_from_Buffers_z(matrix,ibuf,h_buf,matrix_diagonal,    &
                                               non_zero)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:,:)            :: matrix    
  INTEGER, DIMENSION(2,*)               :: ibuf
  COMPLEX*16, DIMENSION(:)              :: h_buf
  COMPLEX*16, DIMENSION(:)              :: matrix_diagonal
  INTEGER                               :: non_zero
  INTEGER                               :: i
  INTEGER                               :: ind_i
  INTEGER                               :: ind_j
  CHARACTER(LEN=80)                     :: title
  matrix(:,:) = 0.d0
  DO i=1,non_zero
     ind_i = ibuf(1,i)
     ind_j = ibuf(2,i)
     matrix(ind_i,ind_j) = h_buf(i)
     matrix(ind_j,ind_i) = conjg(h_buf(i))
  END DO
  DO i = 1,n3d
     matrix(i,i) = matrix_diagonal(i)
  END DO
  IF (print_internal_matrices) THEN
      title= 'Printing Matrix from Fill_Square_Matrix_from_Buffers Subroutine'     
      Call prntcmn(title,matrix,n3d,n3d,n3d,n3d,iout,'e')
  END IF
!***********************************************************************
  END SUBROUTINE Fill_Square_Matrix_from_Buffers_z
!***********************************************************************
!***********************************************************************
!Deck Fill_Buffers_from_Triangular_Matrix_d
!***begin prologue     Fill_Buffers_from_Triangular_Matrix_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill Buffers
!***
!***references
!***routines called
!***end prologue       Fill_Buffers_from_Triangular_Matrix using two index
!***                   non zero method of storage.
!
  SUBROUTINE Fill_Buffers_from_Triangular_Matrix_d(matrix,ibuf,h_buf,matrix_diagonal,    &
                                                   drop,non_zero )
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                  :: matrix
  INTEGER, DIMENSION(2,*)               :: ibuf
  REAL*8, DIMENSION(:)                  :: h_buf
  REAL*8, DIMENSION(:)                  :: matrix_diagonal
  REAL*8                                :: drop
  REAL*8                                :: tmp 
  REAL*8                                :: matrix_el 
  INTEGER                               :: non_zero
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: ind_i
  INTEGER                               :: ind_j
  INTEGER                               :: ij
  INTEGER                               :: count
  CHARACTER(LEN=80)                     :: title
  count = 0
  ij = 0
  DO i=1, n3d
     DO j = 1, i - 1
        ij  = ij + 1
        matrix_el = matrix(ij)
        tmp = abs(matrix_el)
        IF(tmp >= drop) THEN
           count = count + 1
           ibuf(1,count) = i
           ibuf(2,count) = j
           h_buf(count) = matrix_el
        END IF
     END DO
     ij = ij + 1
     matrix_diagonal(i) = matrix_el
  END DO
  non_zero = count
  IF (print_buffers) THEN
     DO i = 1, non_zero
        write(iout,*) 'I = ',ibuf(1,i),'  J = ', ibuf(2,i),'  Val = ', h_buf(i)
     END DO
     write(iout,*) ' Diagonal = ', matrix_diagonal(1:n3d)
  END IF
!***********************************************************************
  END SUBROUTINE Fill_Buffers_from_Triangular_Matrix_d
!***********************************************************************
!***********************************************************************
!Deck Fill_Buffers_from_Triangular_Matrix_z
!***begin prologue     Fill_Buffers_from_Triangular_Matrix_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill Buffers
!***
!***references
!***routines called
!***end prologue       Fill_Buffers_from_Triangular_Matrix_z
!
  SUBROUTINE Fill_Buffers_from_Triangular_Matrix_z(matrix,ibuf,h_buf,matrix_diagonal,          &
                                                   drop,non_zero )
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)              :: matrix
  INTEGER, DIMENSION(2,*)               :: ibuf
  COMPLEX*16, DIMENSION(:)              :: h_buf
  COMPLEX*16, DIMENSION(:)              :: matrix_diagonal
  REAL*8                                :: drop
  REAL*8                                :: tmp
  COMPLEX*16                            :: matrix_el
  INTEGER                               :: non_zero
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: ij
  INTEGER                               :: count
  count = 0
  ij = 0
  DO i=1, n3d
     DO j = 1, i - 1
        ij  = ij + 1
        matrix_el = matrix(ij)
        tmp = abs(matrix_el)
        IF(tmp >= drop) THEN
           count = count + 1
           ibuf(1,count) = i
           ibuf(2,count) = j
           h_buf(count) = matrix_el
        END IF
     END DO
     ij = ij + 1
     matrix_diagonal(i) = matrix_el
  END DO
  non_zero = count
  IF (print_buffers) THEN
     DO i = 1, non_zero
        write(iout,*) 'I = ',ibuf(1,i),'  J = ', ibuf(2,i),'  Val = ', h_buf(i)
     END DO
     write(iout,*) ' Diagonal = ', matrix_diagonal(1:n3d)
  END IF
!***********************************************************************
  END SUBROUTINE Fill_Buffers_from_Triangular_Matrix_z
!***********************************************************************
!***********************************************************************
!Deck Fill_Square_Matrix_from_Disk_d
!***begin prologue     Fill_Square_Matrix_from_Disk_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill_Square
!***
!***references
!***routines called
!***end prologue       Fill_Square_Matrix_from_Disk_f_d
!
  SUBROUTINE Fill_Square_Matrix_From_Disk_d(matrix,non_zero,unit)
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)                :: matrix
  REAL*8                                :: matrix_el
  INTEGER                               :: non_zero
  INTEGER                               :: i
  INTEGER                               :: ind_i
  INTEGER                               :: ind_j
  INTEGER                               :: unit
  CHARACTER(LEN=80)                     :: title 
  matrix(:,:) = 0.d0
  READ (unit) non_zero
  DO i=1,non_zero
     READ(unit) ind_i, ind_j, matrix_el
     matrix(ind_i,ind_j) = matrix_el
     matrix(ind_j,ind_i) = matrix_el
  END DO
  IF (print_internal_matrices) THEN
      title= 'Printing Matrix from Fill_Square_Matrix_from_Buffers Subroutine'     
      Call prntfmn(title,matrix,n3d,n3d,n3d,n3d,iout,'e')
  END IF
!***********************************************************************
  END SUBROUTINE Fill_Square_Matrix_from_Disk_d
!***********************************************************************
!***********************************************************************
!Deck Fill_Square_Matrix_from_Disk_z
!***begin prologue     Fill_Square_Matrix_from_Disk_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill_Square
!***
!***references
!***routines called
!***end prologue       Fill_Square_Matrix_from_Disk_z
!
  SUBROUTINE Fill_Square_Matrix_from_Disk_z(matrix,non_zero,unit)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:,:)            :: matrix
  COMPLEX*16                            :: matrix_el
  INTEGER                               :: non_zero
  INTEGER                               :: i
  INTEGER                               :: ind_i
  INTEGER                               :: ind_j
  INTEGER                               :: unit
  CHARACTER(LEN=80)                     :: title 
  matrix(:,:) = 0.d0
  READ (unit) non_zero
  DO i=1,non_zero
     READ(unit) ind_i, ind_j, matrix_el
     matrix(ind_i,ind_j) = matrix_el
     matrix(ind_j,ind_i) = conjg( matrix_el )
  END DO
  IF (print_internal_matrices) THEN
      title= 'Printing Matrix from Fill_Square_Matrix_from_Buffers Subroutine'     
      Call prntcmn(title,matrix,n3d,n3d,n3d,n3d,iout,'e')
  END IF
!***********************************************************************
  END SUBROUTINE Fill_Square_Matrix_from_Disk_z
!***********************************************************************
!***********************************************************************
!Deck Fill_Triangular_Matrix_from_Disk_d
!***begin prologue     Fill_Triangular_Matrix_from_Disk_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill_Triangular
!***
!***references
!***routines called
!***end prologue       Fill_Triangular_Matrix_from_Disk_f_d
!
  SUBROUTINE Fill_Triangular_Matrix_From_Disk_d(matrix,non_zero,unit)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                  :: matrix
  REAL*8                                :: matrix_el
  INTEGER                               :: non_zero
  INTEGER                               :: i
  INTEGER                               :: ij
  INTEGER                               :: ind_i
  INTEGER                               :: ind_j
  INTEGER                               :: unit
  matrix(:) = 0.d0
  READ (unit) non_zero
  IF(print_input_matrices) THEN
     DO i=1,non_zero
        READ(unit) ind_i, ind_j, matrix_el
        write(iout,*) 'I = ',ind_i,'  J = ', ind_j,'  Val = ', matrix_el
        ij = ind_i * (ind_i - 1 ) / 2 + ind_j
        matrix(ij) = matrix_el
     END DO
  ELSE
     DO i=1,non_zero
        READ(unit) ind_i, ind_j, matrix_el
        ij = ind_i * (ind_i - 1 ) / 2 + ind_j
        matrix(ij) = matrix_el
     END DO
  END IF
!***********************************************************************
  END SUBROUTINE Fill_Triangular_Matrix_from_Disk_d
!***********************************************************************
!***********************************************************************
!Deck Fill_Triangular_Matrix_from_Disk_z
!***begin prologue     Fill_Triangular_Matrix_from_Disk_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill_Triangular
!***
!***references
!***routines called
!***end prologue       Fill_Triangular_Matrix_from_Disk_z
!
  SUBROUTINE Fill_Triangular_Matrix_from_Disk_z(matrix,non_zero,unit)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)              :: matrix
  COMPLEX*16                            :: matrix_el
  INTEGER                               :: non_zero
  INTEGER                               :: i
  INTEGER                               :: ij
  INTEGER                               :: ind_i
  INTEGER                               :: ind_j
  INTEGER                               :: unit
  matrix(:) = 0.d0
  READ (unit) non_zero
  IF(print_input_matrices) THEN
     DO i=1,non_zero
        READ(unit) ind_i, ind_j, matrix_el
        write(iout,*) 'I = ',ind_i,'  J = ', ind_j,'  Val = ', matrix_el
        ij = ind_i * (ind_i - 1 ) / 2 + ind_j
        matrix(ij) = matrix_el
     END DO
  ELSE
     DO i=1,non_zero
        READ(unit) ind_i, ind_j, matrix_el
        ij = ind_i * (ind_i - 1 ) / 2 + ind_j
        matrix(ij) = matrix_el
     END DO
  END IF
!***********************************************************************
  END SUBROUTINE Fill_Triangular_Matrix_from_Disk_z
!***********************************************************************
!***********************************************************************
!Deck Fill_Buffers_from_Disk_d
!***begin prologue     Fill_Buffers_from_Disk_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill Buffers from a disk.
!***
!***references
!***routines called
!***end prologue       Fill_Buffers_from_Disk_d
!
  SUBROUTINE Fill_Buffers_from_Disk_d(ibuf,h_buf,matrix_diagonal,drop,non_zero,unit)
  IMPLICIT NONE
  INTEGER, DIMENSION(2,*)               :: ibuf
  REAL*8, DIMENSION(:)                  :: h_buf
  REAL*8, DIMENSION(:)                  :: matrix_diagonal
  REAL*8                                :: drop
  REAL*8                                :: tmp
  REAL*8                                :: matrix_el
  INTEGER                               :: non_zero
  INTEGER                               :: i
  INTEGER                               :: ind_i
  INTEGER                               :: ind_j
  INTEGER                               :: count
  INTEGER                               :: unit
  count = 0
  DO i=1,non_zero
     READ(unit) ind_i, ind_j, matrix_el
     tmp = abs(matrix_el)
     IF(tmp >= drop) THEN
         IF (ind_i == ind_j) THEN
             matrix_diagonal(ind_i) = matrix_el
         ELSE
             count = count + 1
             ibuf(1,count) = ind_i
             ibuf(2,count) = ind_j
             h_buf(count) = matrix_el
         END IF
      END IF
   END DO
  non_zero = count
  IF (print_buffers) THEN
     DO i = 1, non_zero
        write(iout,*) 'I = ',ibuf(1,i),'  J = ', ibuf(2,i),'  Val = ', h_buf(i)
     END DO
     write(iout,*) ' Diagonal = ', matrix_diagonal(1:n3d)
  END IF
!***********************************************************************
  END SUBROUTINE Fill_Buffers_from_Disk_d
!***********************************************************************
!***********************************************************************
!Deck Fill_Buffers_from_Disk_z
!***begin prologue     Fill_Buffers_from_Disk_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill Buffers
!***
!***references
!***routines called
!***end prologue       Fill_Buffers_from_Disk_z
!
  SUBROUTINE Fill_Buffers_from_Disk_z(ibuf,h_buf,matrix_diagonal,drop,non_zero,unit)
  IMPLICIT NONE
  INTEGER, DIMENSION(2,*)               :: ibuf
  COMPLEX*16, DIMENSION(:)              :: h_buf
  COMPLEX*16, DIMENSION(:)              :: matrix_diagonal
  REAL*8                                :: drop
  REAL*8                                :: tmp
  COMPLEX*16                            :: matrix_el
  INTEGER                               :: non_zero
  INTEGER                               :: i
  INTEGER                               :: ind_i
  INTEGER                               :: ind_j
  INTEGER                               :: count
  INTEGER                               :: unit
  count = 0
  DO i=1,non_zero
     READ(unit) ind_i, ind_j, matrix_el
     tmp = abs(matrix_el)
     IF(tmp >= drop) THEN
         IF (ind_i == ind_j) THEN
             matrix_diagonal(ind_i) = matrix_el
         ELSE
             count = count + 1
             ibuf(1,count) = ind_i
             ibuf(2,count) = ind_j
             h_buf(count) = matrix_el
         END IF
      END IF
   END DO
  non_zero = count
  IF (print_buffers) THEN
     DO i = 1, non_zero
        write(iout,*) 'I = ',ibuf(1,i),'  J = ', ibuf(2,i),'  Val = ', h_buf(i)
     END DO
     write(iout,*) ' Diagonal = ', matrix_diagonal(1:n3d)
  END IF
!***********************************************************************
  END SUBROUTINE Fill_Buffers_from_Disk_z
!***********************************************************************
!***********************************************************************
!Deck Fill_Square_Matrix_from_Input_d
!***begin prologue     Fill_Square_Matrix_from_Input_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill_Square
!***
!***references
!***routines called
!***end prologue       Fill_Square_Matrix_from_Input_d
!
  SUBROUTINE Fill_Square_Matrix_From_Input_d(keyword,matrix,n)
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)                :: matrix
  INTEGER                               :: n
  INTEGER                               :: i
  LOGICAL                               :: posinp
  CHARACTER(LEN=80)                     :: cpass 
  CHARACTER(LEN=*)                      :: keyword 
  IF ( posinp(keyword,cpass ) ) THEN
       Write(iout,*)
       Write(iout,*) ' Reading in Square Matrix'
       Write(iout,*) ' Keyword = ',keyword
       Write(iout,*)
       DO i=1,n
          READ(inp,*) matrix(i,:)
          WRITE(iout,1) i , matrix(i,:)
       END DO
  END IF
  Write(iout,*)
1 Format(/,5x,'Row = ',i5,(/,5x,5f15.8))
!***********************************************************************
  END SUBROUTINE Fill_Square_Matrix_from_Input_d
!***********************************************************************
!***********************************************************************
!Deck Fill_Square_Matrix_from_Input_z
!***begin prologue     Fill_Square_Matrix_from_Input_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill_Square
!***
!***references
!***routines called
!***end prologue       Fill_Square_Matrix_from_Input_z
!
  SUBROUTINE Fill_Square_Matrix_from_Input_z(keyword,matrix,n)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:,:)            :: matrix
  INTEGER                               :: n
  INTEGER                               :: i
  CHARACTER(LEN=*)                      :: keyword 
  CHARACTER(LEN=80)                     :: cpass 
  LOGICAL                               :: posinp
  IF ( posinp(keyword,cpass ) ) THEN
       Write(iout,*)
       Write(iout,*) ' Reading in Square Matrix'
       Write(iout,*) ' Keyword = ',keyword
       Write(iout,*)
       DO i=1,n
          READ(inp,*) matrix(i,:)
          WRITE(iout,1) i , matrix(i,:)
       END DO
  END IF
  Write(iout,*)
1 Format(/,5x,'Row = ',i5,(/,5x,5f15.8))
!***********************************************************************
  END SUBROUTINE Fill_Square_Matrix_from_Input_z
!***********************************************************************
!***********************************************************************
!Deck Fill_Triangular_Matrix_from_Input_d
!***begin prologue     Fill_Triangular_Matrix_from_Input_d
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill_Triangular
!***
!***references
!***routines called
!***end prologue       Fill_Triangular_Matrix_from_Input_f_d
!
  SUBROUTINE Fill_Triangular_Matrix_From_Input_d(keyword,matrix,n)
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                  :: matrix
  INTEGER                               :: n
  INTEGER                               :: i
  INTEGER                               :: begin
  INTEGER                               :: end
  LOGICAL                               :: posinp
  CHARACTER(LEN=*)                      :: keyword 
  CHARACTER(LEN=80)                     :: cpass 
  IF ( posinp(keyword,cpass ) ) THEN
       Write(iout,*)
       Write(iout,*) ' Reading in Triangular Matrix'
       Write(iout,*) ' Keyword = ',keyword
       Write(iout,*)
       begin = 0
       end = 0
       DO i=1,n
          begin = end + 1
          end = end + i 
          READ(inp,*) matrix(begin:end)
          WRITE(iout,1) i, matrix(begin:end)
       END DO
  END IF
  Write(iout,*)
1 Format(/,5x,'Row = ',i5,(/,5x,5f15.8))
!***********************************************************************
  END SUBROUTINE Fill_Triangular_Matrix_from_Input_d
!***********************************************************************
!***********************************************************************
!Deck Fill_Triangular_Matrix_from_Input_z
!***begin prologue     Fill_Triangular_Matrix_from_Input_z
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill_Triangular
!***
!***references
!***routines called
!***end prologue       Fill_Triangular_Matrix_from_Input_z
!
 SUBROUTINE Fill_Triangular_Matrix_from_Input_z(keyword,matrix,n)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)              :: matrix
  INTEGER                               :: n
  INTEGER                               :: i
  INTEGER                               :: begin
  INTEGER                               :: end
  LOGICAL                               :: posinp
  CHARACTER(LEN=80)                     :: cpass 
  CHARACTER(LEN=*)                      :: keyword 
  IF ( posinp(keyword,cpass ) ) THEN
       Write(iout,*)
       Write(iout,*) ' Reading in Triangular Matrix'
       Write(iout,*) ' Keyword = ',keyword
       Write(iout,*)
       begin = 0
       end = 0
       DO i=1,n
          begin = end + 1
          end = end + i 
          READ(inp,*) matrix(begin:end)
          WRITE(iout,1) i, matrix(begin:end)
       END DO
  END IF
  Write(iout,*)
1 Format(/,5x,'Row = ',i5,(/,5x,5f15.8))
!***********************************************************************
  END SUBROUTINE Fill_Triangular_Matrix_from_Input_z
!***********************************************************************
  END  MODULE Matrix_Utility_Module
!***********************************************************************
!***********************************************************************
