!***********************************************************************
                           MODULE Pack_Hamiltonian_Module
                            USE input_output
                            USE Pack_Global,                                             &
                                ONLY : drop_tol
                            USE Global_Time_Propagation_Module,                          &
                                ONLY : print_packed_matrices
                            IMPLICIT NONE
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          INTERFACE h_pack_upper
                   MODULE PROCEDURE h_pack_upper_d,                                      &
                                    h_pack_upper_z
                          END INTERFACE h_pack_upper
!
                          INTERFACE h_pack_lower
                   MODULE PROCEDURE h_pack_lower_d,                                      &
                                    h_pack_lower_z
                          END INTERFACE h_pack_lower
!
                          INTERFACE h_pack_random
                   MODULE PROCEDURE h_pack_random_d,                                     &
                                    h_pack_random_z
                          END INTERFACE h_pack_random
!
                          INTERFACE pack_matrix
                   MODULE PROCEDURE pack_matrix_d,                                       &
                                    pack_matrix_z
                          END INTERFACE pack_matrix
!
                          INTERFACE Pack_Rectangular
                   MODULE PROCEDURE Pack_Rectangular_d,                                  &
                                    Pack_Rectangular_z
                          END INTERFACE Pack_Rectangular
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           Contains
!***********************************************************************
!***********************************************************************
!***begin prologue     h_pack_upper_d
!***date written       000710   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           one-dim
!***author             schneider, barry (nsf)
!***source             
!***purpose            Pack non zero hamiltonian matrix elements and indices. 
!***                   What is stored are the non zero elements of each column
!***                   of the upper triangle.  The array packed_columns has
!***                   these.  The row index is stored is row_index and the
!***                   number of non zero elements in the column in non_zero_column.
!***                   Storage is a11, a12, a22, a13, a23, a33 etc.
!***references         
!
!***routines called    
!***end prologue       h_pack_upper_d
  Subroutine h_pack_upper_d(upper,matrix_diagonal,packed_columns,               &
                            non_zero_columns,row_index,number,n,type)
  IMPLICIT NONE
  REAL*8,   DIMENSION(:)       :: upper
  REAL*8,   DIMENSION(:)       :: matrix_diagonal
  REAL*8,   DIMENSION(:)       :: packed_columns
  INTEGER,  DIMENSION(:)       :: non_zero_columns
  INTEGER,  DIMENSION(:)       :: row_index
  INTEGER                      :: n
  INTEGER                      :: col_ind
  INTEGER                      :: row_ind
  INTEGER                      :: num
  INTEGER                      :: number
  INTEGER                      :: count
  CHARACTER(LEN=*)             :: type
  CHARACTER(LEN=80)            :: title
!
!     pack all non-zero elements
! 
  count = 0  
  non_zero_columns(:) = 0
  num = 0
  write(iout,1) type
  DO col_ind = 1 , n
     DO row_ind = 1 , col_ind-1
        num = num + 1
        IF(abs(upper(num)) >= drop_tol) THEN 
           non_zero_columns(col_ind) = non_zero_columns(col_ind) + 1
           count = count + 1
           row_index(count) = row_ind
           packed_columns(count) = upper(num)
        END IF
     END DO
     num = num + 1
     matrix_diagonal(col_ind) = upper(num)
  END DO
  number= count
!
!
  IF(print_packed_matrices) THEN
     WRITE(iout,2) n*(n+1)/2, number, drop_tol
     title='Matrix Diagonal'
     call prntfmn(title,matrix_diagonal,n,1,n,1,iout,'e')
     title='non-zero columns'  
     call prntim(title,non_zero_columns,n,1,n,1,iout)
     title='row_index'  
     call prntim(title,row_index,count,1,count,1,iout)
     title='packed_columns'  
     call prntfmn(title,packed_columns,count,1,count,1,iout,'e')
  END IF
1 FORMAT(/,10x,'Packed Matrix = ', a16,'  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Possible Non-Zero Elements = ', i15,                              &
           10x,'Actual Non-Zero Elements   = ',i15,                               &
         /10x, 'Drop_Tolerance             = ',e15.8)
END SUBROUTINE h_pack_upper_d
!***********************************************************************
!***********************************************************************
!***begin prologue     h_pack_upper_z
!***date written       000710   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           one-dim
!***author             schneider, barry (nsf)
!***source             
!***purpose            pack non zero hamiltonian matrix elements and indices 
!***                   
!***references         
!
!***routines called    
!***end prologue       h_pack_upper_z
 Subroutine h_pack_upper_z(upper,matrix_diagonal,packed_columns,               &
                           non_zero_columns,row_index,number,n,type)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)     :: upper
  COMPLEX*16, DIMENSION(:)     :: matrix_diagonal
  COMPLEX*16, DIMENSION(:)     :: packed_columns
  INTEGER,  DIMENSION(:)       :: non_zero_columns
  INTEGER,  DIMENSION(:)       :: row_index
  INTEGER                      :: n
  INTEGER                      :: col_ind
  INTEGER                      :: row_ind
  INTEGER                      :: num
  INTEGER                      :: number
  INTEGER                      :: count
  CHARACTER(LEN=80)            :: title
  CHARACTER(LEN=*)             :: type
!
!     pack all non-zero elements
! 
  count=0
  non_zero_columns(:) = 0
  num = 0
  write(iout,1) type
  DO col_ind = 1 , n
     DO row_ind = 1 , col_ind-1
        num = num + 1
        IF(abs(upper(num)) >= drop_tol) THEN 
           non_zero_columns(col_ind) = non_zero_columns(col_ind) + 1
           count = count + 1
           row_index(count) = row_ind
           packed_columns(count) = upper(num)
        END IF
     END DO
     num = num + 1
     matrix_diagonal(col_ind) = upper(num)
  END DO
  number = count
! 
! 
  IF(print_packed_matrices) THEN
     WRITE(iout,2) n*(n+1)/2, number, drop_tol
     title='Matrix Diagonal'
     call prntcmn(title,matrix_diagonal,n,1,n,1,iout,'e')
     title='non-zero columns'  
     call prntim(title,non_zero_columns,n,1,n,1,iout)
     title='row_index'  
     call prntim(title,row_index,count,1,count,1,iout)
     title='packed_columns'  
     call prntcmn(title,packed_columns,count,1,count,1,iout,'e')
  END IF
1 FORMAT(/,10x,'Packed Matrix = ', a16,'  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Possible Non-Zero Elements = ', i15,                              &
           10x,'Actual Non-Zero Elements   = ',i15,                               &
         /10x, 'Drop_Tolerance             = ',e15.8)
END SUBROUTINE h_pack_upper_z
!***********************************************************************
!***********************************************************************
!***begin prologue     h_pack_lower_d
!***date written       000710   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           one-dim
!***author             schneider, barry (nsf)
!***source             
!***purpose            Pack non zero hamiltonian matrix elements and indices 
!***                   The same idea is used as in the upper routine but here
!***                   applied to the lower triangle.
!***                   Storage is a11, a21, a31, .... an1, a22, a32, a33 etc.
!***references         
!
!***routines called    
!***end prologue       h_pack_lower_d
  Subroutine h_pack_lower_d(lower,matrix_diagonal,packed_columns,               &
                            non_zero_columns,row_index,number,n,type)
  IMPLICIT NONE
  REAL*8,   DIMENSION(:)       :: lower
  REAL*8,   DIMENSION(:)       :: matrix_diagonal
  REAL*8,   DIMENSION(:)       :: packed_columns
  INTEGER,  DIMENSION(:)       :: non_zero_columns
  INTEGER,  DIMENSION(:)       :: row_index
  INTEGER                      :: n
  INTEGER                      :: row_ind
  INTEGER                      :: col_ind
  INTEGER                      :: num
  INTEGER                      :: number
  INTEGER                      :: count
  CHARACTER(LEN=80)            :: title
  CHARACTER(LEN=*)             :: type
!
!     pack all non-zero elements
! 
  count = 0  
  non_zero_columns(:) = 0
  num = 0
  write(iout,1) type
  DO col_ind = 1 , n
     num = num + 1
     DO row_ind = col_ind + 1, n
        num = num + 1
        IF(abs(lower(num)) >= drop_tol) THEN 
           non_zero_columns(col_ind) = non_zero_columns(col_ind) + 1
           count = count + 1
           row_index(count) = row_ind
           packed_columns(count) = lower(num)
        END IF
     END DO
  END DO
  number= count
!
!
  IF(print_packed_matrices) THEN
     WRITE(iout,2) n*(n+1)/2, number, drop_tol
     title='Matrix Diagonal'
     call prntfmn(title,matrix_diagonal,n,1,n,1,iout,'e')
     title='non-zero columns'  
     call prntim(title,non_zero_columns,n,1,n,1,iout)
     title='row_index'  
     call prntim(title,row_index,count,1,count,1,iout)
     title='packed_columns'  
     call prntfmn(title,packed_columns,count,1,count,1,iout,'e')
  END IF
1 FORMAT(/,10x,'Packed Matrix = ', a16,' All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Possible Non-Zero Elements = ', i15,                              &
           10x,'Actual Non-Zero Elements   = ',i15,                               &
         /10x, 'Drop_Tolerance             = ',e15.8)
END SUBROUTINE h_pack_lower_d
!***********************************************************************
!***********************************************************************
!***begin prologue     h_pack_lower_z
!***date written       000710   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           one-dim
!***author             schneider, barry (nsf)
!***source             
!***purpose            pack non zero hamiltonian matrix elements and indices 
!***                   
!***references         
!
!***routines called    
!***end prologue       h_pack_lower_z
  Subroutine h_pack_lower_z(lower,matrix_diagonal,packed_columns,               &
                            non_zero_columns,row_index,number,n,type)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)     :: lower
  COMPLEX*16, DIMENSION(:)     :: matrix_diagonal
  COMPLEX*16, DIMENSION(:)     :: packed_columns
  INTEGER,  DIMENSION(:)       :: non_zero_columns
  INTEGER,  DIMENSION(:)       :: row_index
  INTEGER                      :: n
  INTEGER                      :: row_ind
  INTEGER                      :: col_ind
  INTEGER                      :: num
  INTEGER                      :: number
  INTEGER                      :: count
  CHARACTER(LEN=80)            :: title
  CHARACTER(LEN=*)             :: type
!
!     pack all non-zero elements
! 
  count=0
  non_zero_columns(:) = 0
  num = 0
  write(iout,1) type
  DO col_ind = 1 , n
     num = num + 1
     DO row_ind = col_ind + 1 , n
        num = num + 1
        IF(abs(lower(num)) >= drop_tol) THEN 
           non_zero_columns(col_ind) = non_zero_columns(col_ind) + 1
           count = count + 1
           row_index(count) = row_ind
           packed_columns(count) = lower(num)
        END IF
     END DO
  END DO
  number = count
! 
! 
  IF(print_packed_matrices) THEN
     WRITE(iout,2) n*(n+1)/2, number, drop_tol
     title='Matrix Diagonal'
     call prntcmn(title,matrix_diagonal,n,1,n,1,iout,'e')
     title='non-zero columns'  
     call prntim(title,non_zero_columns,n,1,n,1,iout)
     title='row_index'  
     call prntim(title,row_index,count,1,count,1,iout)
     title='packed_columns'  
     call prntcmn(title,packed_columns,count,1,count,1,iout,'e')
  END IF
1 FORMAT(/,10x,'Packed Matrix = ', a16,'  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Possible Non-Zero Elements = ', i15,                              &
           10x,'Actual Non-Zero Elements   = ',i15,                               &
         /10x, 'Drop_Tolerance             = ',e15.8)
END SUBROUTINE h_pack_lower_z
!***********************************************************************
!***********************************************************************
!*deck  @(#)h_pack_random_d
  subroutine h_pack_random_d(ham,matrix_diagonal,hbuf,ibuf,number,n,type)
!***revision date      910606   (yymmdd)
!***end prologue       h_pack_random_d
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)            :: ham
  REAL*8, DIMENSION(:)              :: matrix_diagonal
  REAL*8, DIMENSION(:)              :: hbuf
  INTEGER, DIMENSION(:,:)           :: ibuf
  INTEGER                           :: n
  INTEGER                           :: i
  INTEGER                           :: j
  INTEGER                           :: number
  INTEGER                           :: ntot
  CHARACTER(LEN=*)                  :: type
  number=0
  DO i=1,n
     DO j=1, i-1
        IF(abs(ham(i,j)).gt.drop_tol) then
           number=number+1
           ibuf(number,1) = i
           ibuf(number,2) = j
           hbuf(number)=ham(i,j)
        END IF
     END DO
  END DO
  DO i=1,n
     matrix_diagonal(i)=ham(i,i)
  END DO
  ntot = n*(n+1)/2
  WRITE(iout,1) ntot, number, drop_tol
1 FORMAT(/,10x,'Possible Non-Zero Elements =', i12,                              &
         /10x,'Actual Non-Zero Elements    =',i12,                               &
         /10x,'Drop_Tolerance              =',d15.8)
  END SUBROUTINE h_pack_random_d
!***********************************************************************
!***********************************************************************
!*deck  @(#)h_pack_random_z
  subroutine h_pack_random_z(ham,matrix_diagonal,hbuf,ibuf,n,type)
!***revision date      910606   (yymmdd)
!***end prologue       h_pack_random_z
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:,:)            :: ham
  COMPLEX*16, DIMENSION(:)              :: matrix_diagonal
  COMPLEX*16, DIMENSION(:)              :: hbuf
  INTEGER,    DIMENSION(:,:)            :: ibuf
  INTEGER                               :: n
  INTEGER                               :: i
  INTEGER                               :: j
  INTEGER                               :: number
  INTEGER                               :: ntot
  CHARACTER(LEN=*)                      :: type
  write(iout,1) type
  number=0
  DO i=1,n
     DO j=1, i-1
        IF(abs(ham(i,j)).gt.drop_tol) then
           number=number+1
           ibuf(number,1) = i
           ibuf(number,2) = j
           hbuf(number)=ham(i,j)
        END IF
     END DO
  END DO
  DO i=1,n
     matrix_diagonal(i)=ham(i,i)
  END DO
  ntot = n*(n+1)/2
  WRITE(iout,2) ntot, drop_tol, number               
1 FORMAT(/,10x,'Packed Matrix = ', a16,'  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Possible Non Zero Matrix Elements = ',i12,                         &
           10x,'Actual Non Zero Matrix Elements   = ',i12,                         &
         /10x, 'Drop_Tolerance                    = ',e15.8)
  END SUBROUTINE h_pack_random_z
!***********************************************************************
!***********************************************************************
!***begin prologue     pack_rectangular_d
!***date written       000710   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           one-dim
!***author             schneider, barry (nsf)
!***source             
!***purpose            Pack non zero hamiltonian matrix elements and indices. 
!***                   What is stored are the non zero elements of each column
!***                   of the upper triangle.  The array packed_columns has
!***                   these.  The row index is stored in row_index and the
!***                   number of non zero elements in the column in non_zero_column.
!***                   Storage is a11, a12, a22, a13, a23, a33 etc.
!***references         
!
!***routines called    
!***end prologue       pack_rectangular_d
  Subroutine pack_rectangular_d(matrix,packed_columns,                           &
                                non_zero_columns,row_index,number,               &
                                n_row,n_col,type)
  IMPLICIT NONE
  REAL*8,   DIMENSION(:,:)     :: matrix
  REAL*8,   DIMENSION(:)       :: packed_columns
  INTEGER,  DIMENSION(:)       :: non_zero_columns
  INTEGER,  DIMENSION(:)       :: row_index
  INTEGER                      :: n_row
  INTEGER                      :: n_col
  INTEGER                      :: i
  INTEGER                      :: j
  INTEGER                      :: number
  INTEGER                      :: count
  CHARACTER(LEN=*)             :: type
  CHARACTER(LEN=80)            :: title
!
!     pack all non-zero elements
! 
  count = 0  
  non_zero_columns(:) = 0
  write(iout,1) type
  DO j=1,n_col
     DO i=1,n_row
        IF(abs(matrix(i,j)) >= drop_tol) THEN 
           non_zero_columns(j) = non_zero_columns(j) + 1
           count = count + 1
           row_index(count) = i
           packed_columns(count) = matrix(i,j)
        END IF
     END DO
  END DO
  number= count
!
!
!
  IF(print_packed_matrices) THEN
     WRITE(iout,2) n_col*n_row, number, drop_tol
     title='non-zero columns'  
     call prntim(title,non_zero_columns,n_col,1,n_col,1,iout)
     title='row_index'  
     call prntim(title,row_index,count,1,count,1,iout)
     title='packed_columns'  
     call prntfmn(title,packed_columns,count,1,count,1,iout,'e')
  END IF
1 FORMAT(/,10x,'Packed Matrix = ', a16,'  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Possible Non-Zero Elements = ', i15,                              &
           10x,'Actual Non-Zero Elements   = ',i15,                               &
         /10x, 'Drop_Tolerance             = ',e15.8)
END SUBROUTINE pack_rectangular_d
!***********************************************************************
!***********************************************************************
!***begin prologue     pack_rectangular_z
!***date written       000710   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           one-dim
!***author             schneider, barry (nsf)
!***source             
!***purpose            pack non zero hamiltonian matrix elements and indices 
!***                   
!***references         
!
!***routines called    
!***end prologue       pack_rectangular_z
 Subroutine pack_rectangular_z(matrix,packed_columns,                            &
                               non_zero_columns,row_index,number,                &
                               n_row,n_col,type)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:,:)   :: matrix
  COMPLEX*16, DIMENSION(:)     :: packed_columns
  INTEGER,  DIMENSION(:)       :: non_zero_columns
  INTEGER,  DIMENSION(:)       :: row_index
  INTEGER                      :: n_row
  INTEGER                      :: n_col
  INTEGER                      :: i
  INTEGER                      :: j
  INTEGER                      :: number
  INTEGER                      :: count
  CHARACTER(LEN=80)            :: title
  CHARACTER(LEN=*)             :: type
!
!     pack all non-zero elements
! 
  count = 0  
  non_zero_columns(:) = 0
  write(iout,1) type
  DO j=1,n_col
     DO i=1,n_row
        IF(abs(matrix(i,j)) >= drop_tol) THEN 
           non_zero_columns(j) = non_zero_columns(j) + 1
           count = count + 1
           row_index(count) = i
           packed_columns(count) = matrix(i,j)
        END IF
     END DO
  END DO
  number = count
! 
! 
  IF(print_packed_matrices) THEN
     WRITE(iout,2) n_col*n_row, number, drop_tol
     title='non-zero columns'  
     call prntim(title,non_zero_columns,n_col,1,n_col,1,iout)
     title='row_index'  
     call prntim(title,row_index,count,1,count,1,iout)
     title='packed_columns'  
     call prntcmn(title,packed_columns,count,1,count,1,iout,'e')
  END IF
1 FORMAT(/,10x,'Packed Matrix = ', a16,'  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Possible Non-Zero Elements = ', i15,                              &
           10x,'Actual Non-Zero Elements   = ',i15,                               &
         /10x, 'Drop_Tolerance             = ',e15.8)
END SUBROUTINE pack_rectangular_z
!***********************************************************************
!***********************************************************************
  Subroutine pack_matrix_d(matrix,matrix_diagonal,packed_columns,                 &
                           non_zero_columns,row_index,number,n,type)

!***revision date      910606   (yymmdd)
!***end prologue       pack_matrix_d
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)            :: matrix
  REAL*8, DIMENSION(:)              :: matrix_diagonal
  REAL*8, DIMENSION(:)              :: packed_columns
  REAL*8                            :: mat_el
  INTEGER, DIMENSION(:)             :: non_zero_columns
  INTEGER, DIMENSION(:)             :: row_index
  INTEGER                           :: n
  INTEGER                           :: i
  INTEGER                           :: j
  INTEGER                           :: number
  INTEGER                           :: num
  CHARACTER(LEN=*)                  :: type
  CHARACTER(LEN=80)                 :: title
  write(iout,1) type
  number=0
  non_zero_columns(:) = 0
  DO j=1,n
     DO i=j+1, n
        mat_el = matrix(i,j)
        IF(abs(mat_el) >= drop_tol) THEN 
           non_zero_columns(j) = non_zero_columns(j) + 1
           number = number + 1
           row_index(number) = i
           packed_columns(number) = mat_el
        END IF
     END DO
  END DO
  DO i=1,n
     matrix_diagonal(i) = matrix(i,i)
  END DO
  num = n*(n+1)/2
  IF(print_packed_matrices) THEN
     WRITE(iout,2) num, number, drop_tol
     title='Matrix Diagonal'
     call prntfmn(title,matrix_diagonal,n,1,n,1,'e')
     title='non-zero columns'  
     call prntim(title,non_zero_columns,n,1,n,1)
     title='row_index'  
     call prntim(title,row_index,number,1,number,1)
     title='packed_columns'  
     call prntfmn(title,packed_columns,number,1,number,1,'e')
  END IF
1 FORMAT(/,10x,'Packed Matrix = ', a16,'  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Possible Non-Zero Elements = ', i12,                              &
           10x,'Actual Non-Zero Elements   = ',i12,                               &
         /10x, 'Drop_Tolerance             = ',e15.8)
  END SUBROUTINE pack_matrix_d
!***********************************************************************
!***********************************************************************
  subroutine pack_matrix_z(matrix,matrix_diagonal,packed_columns,                 &
                           non_zero_columns,row_index,number,n,type)
!***revision date      910606   (yymmdd)
!***end prologue      pack_matrix_z
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:,:)        :: matrix
  COMPLEX*16, DIMENSION(:)          :: matrix_diagonal
  COMPLEX*16, DIMENSION(:)          :: packed_columns
  COMPLEX*16                        :: mat_el
  INTEGER, DIMENSION(:)             :: non_zero_columns
  INTEGER, DIMENSION(:)             :: row_index
  INTEGER                           :: n
  INTEGER                           :: i
  INTEGER                           :: j
  INTEGER                           :: number
  INTEGER                           :: num
  CHARACTER(LEN=*)                  :: type
  CHARACTER(LEN=80)                 :: title
  number=0
  non_zero_columns(:) = 0
  write(iout,1) type
  DO j=1,n
     DO i=j+1, n
        mat_el = matrix(i,j)
        IF(abs(mat_el) >= drop_tol) THEN 
           non_zero_columns(j) = non_zero_columns(j) + 1
           number = number + 1
           row_index(number) = i
           packed_columns(number) = mat_el
        END IF
     END DO
  END DO
  DO i=1,n
     matrix_diagonal(i) = matrix(i,i)
  END DO
  num=n*(n+1)/2
  IF(print_packed_matrices) THEN
     WRITE(iout,2) num, number, drop_tol
     title='Matrix Diagonal'
     call prntcmn(title,matrix_diagonal,n,1,n,1,iout,'e')
     title='non-zero columns'  
     call prntim(title,non_zero_columns,n,1,n,1,iout)
     title='row_index'  
     call prntim(title,row_index,number,1,number,1,iout)
     title='packed_columns'  
     call prntcmn(title,packed_columns,number,1,number,1,iout,'e')
  END IF
1 FORMAT(/,10x,'Packed Matrix = ', a16,'  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Possible Non-Zero Elements = ', i12,                             &
           10x,'Actual Non-Zero Elements   = ',i12,                              &
         /10x, 'Drop_Tolerance             = ',e15.8)
  END SUBROUTINE pack_matrix_z
!***********************************************************************
!***********************************************************************
  END  MODULE Pack_Hamiltonian_Module
!***********************************************************************
!***********************************************************************
