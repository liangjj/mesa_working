!***********************************************************************
                           MODULE Pack_Matrix_Module
                            USE input_output
                            USE dvrprop_global
                            USE Iterative_Global
                            USE Pack_Global
                            IMPLICIT NONE
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          INTERFACE h_pack_upper
                   MODULE PROCEDURE h_pack_upper_d,                             &
                                    h_pack_upper_z
                          END INTERFACE h_pack_upper
!
                          INTERFACE h_pack_lower
                   MODULE PROCEDURE h_pack_lower_d,                             &
                                    h_pack_lower_z
                          END INTERFACE h_pack_lower
!
                          INTERFACE pack_matrix
                   MODULE PROCEDURE pack_matrix_d,                              &
                                    pack_matrix_z
                          END INTERFACE pack_matrix

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
  Subroutine h_pack_upper_d(upper,non_zero_columns,row_index,packed_columns,    &
                            matrix_diagonal,number,prn)
  IMPLICIT NONE
  REAL*8,   DIMENSION(:)       :: upper
  REAL*8,   DIMENSION(:)       :: packed_columns
  REAL*8,   DIMENSION(:)       :: matrix_diagonal
  INTEGER,  DIMENSION(:)       :: non_zero_columns
  INTEGER,  DIMENSION(:)       :: row_index
  INTEGER                      :: i
  INTEGER                      :: j
  INTEGER                      :: num
  INTEGER                      :: number
  INTEGER                      :: count
  LOGICAL                      :: prn
  CHARACTER(LEN=80)            :: title
!
!     pack all non-zero elements
! 
  count = 0  
  non_zero_columns(:) = 0
  num = 0
  write(iout,1)
  DO j=1,n3d
     DO i=1,j-1
        num = num + 1
        IF(abs(upper(num)) >= drop_tol) THEN 
           non_zero_columns(j) = non_zero_columns(j) + 1
           count = count + 1
           row_index(count) = i
           packed_columns(count) = upper(num)
        END IF
     END DO
     num = num + 1
     matrix_diagonal(j) = upper(num)
  END DO
  number= count
!  write(iout,*) number
!  write(iout,*) (non_zero_columns(j), j=1,n3d)
!  write(iout,*) (row_index(j), j=1,count)
!  write(iout,*) (packed_columns(j), j=1,count)
!  write(iout,*) (matrix_diagonal(j), j=1,n3d)
!
!
  WRITE(iout,2) n3d*(n3d+1)/2, number
1 FORMAT(/,10x,'Packed Matrix.  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Possible Non-Zero Elements =', i6,5x, 'Actual Non-Zero Elements =',i6)
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
 Subroutine h_pack_upper_z(upper,non_zero_columns,row_index,packed_columns,    &
                           matrix_diagonal,number,prn)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)     :: upper
  COMPLEX*16, DIMENSION(:)     :: packed_columns
  COMPLEX*16, DIMENSION(:)     :: matrix_diagonal
  INTEGER,  DIMENSION(:)       :: non_zero_columns
  INTEGER,  DIMENSION(:)       :: row_index
  INTEGER                      :: i
  INTEGER                      :: j
  INTEGER                      :: num
  INTEGER                      :: number
  INTEGER                      :: count
  LOGICAL                      :: prn
  CHARACTER(LEN=80)            :: title
!
!     pack all non-zero elements
! 
  count=0
  non_zero_columns(:) = 0
  num = 0
  write(iout,1)
  DO j=1,n3d
     DO i=1,j-1
        num = num + 1
        IF(abs(upper(num)) >= drop_tol) THEN 
           non_zero_columns(j) = non_zero_columns(j) + 1
           count = count + 1
           row_index(count) = i
           packed_columns(count) = upper(num)
        END IF
     END DO
     num = num + 1
     matrix_diagonal(j) = upper(num)
  END DO
  number = count
!  write(iout,*) number
!  write(iout,*) (non_zero_columns(j), j=1,n3d)
!  write(iout,*) (row_index(j), j=1,count)
!  write(iout,*) (packed_columns(j), j=1,count)
!  write(iout,*) (matrix_diagonal(j), j=1,n3d)
! 
! 
  WRITE(iout,2) n3d*(n3d+1)/2, number
1 FORMAT(/,10x,'Packed Matrix.  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Possible Non-Zero Elements =', i6,5x, 'Actual Non-Zero Elements =',i6)
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
  Subroutine h_pack_lower_d(lower,non_zero_columns,row_index,packed_columns,    &
                            matrix_diagonal,number,prn)
  IMPLICIT NONE
  REAL*8,   DIMENSION(:)       :: lower
  REAL*8,   DIMENSION(:)       :: packed_columns
  REAL*8,   DIMENSION(:)       :: matrix_diagonal
  INTEGER,  DIMENSION(:)       :: non_zero_columns
  INTEGER,  DIMENSION(:)       :: row_index
  INTEGER                      :: i
  INTEGER                      :: j
  INTEGER                      :: num
  INTEGER                      :: number
  INTEGER                      :: count
  LOGICAL                      :: prn
  CHARACTER(LEN=80)            :: title
!
!     pack all non-zero elements
! 
  count = 0  
  non_zero_columns(:) = 0
  num = 0
  write(iout,1)
  DO j=1,n3d
     num = num + 1
     DO i=j+1, n3d
        num = num + 1
        IF(abs(lower(num)) >= drop_tol) THEN 
           non_zero_columns(j) = non_zero_columns(j) + 1
           count = count + 1
           row_index(count) = i
           packed_columns(count) = lower(num)
        END IF
     END DO
  END DO
  number= count
!  write(iout,*) number
!  write(iout,*) (non_zero_columns(j), j=1,n3d)
!  write(iout,*) (row_index(j), j=1,count)
!  write(iout,*) (packed_columns(j), j=1,count)
!
!
  WRITE(iout,2) n3d*(n3d+1)/2, number
1 FORMAT(/,10x,'Packed Matrix.  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Possible Non-Zero Elements =', i6,5x, 'Actual Non-Zero Elements =',i6)
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
 Subroutine h_pack_lower_z(lower,non_zero_columns,row_index,packed_columns,    &
                           matrix_diagonal,number,prn)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:)     :: lower
  COMPLEX*16, DIMENSION(:)     :: packed_columns
  COMPLEX*16, DIMENSION(:)     :: matrix_diagonal
  INTEGER,  DIMENSION(:)       :: non_zero_columns
  INTEGER,  DIMENSION(:)       :: row_index
  INTEGER                      :: i
  INTEGER                      :: j
  INTEGER                      :: num
  INTEGER                      :: number
  INTEGER                      :: count
  LOGICAL                      :: prn
  CHARACTER(LEN=80)            :: title
!
!     pack all non-zero elements
! 
  count=0
  non_zero_columns(:) = 0
  num = 0
  write(iout,1)
  DO j=1,n3d
     num = num + 1
     DO i=j+1, n3d
        num = num + 1
        IF(abs(lower(num)) >= drop_tol) THEN 
           non_zero_columns(j) = non_zero_columns(j) + 1
           count = count + 1
           row_index(count) = i
           packed_columns(count) = lower(num)
        END IF
     END DO
  END DO
!  number = count
!  write(iout,*) number
!  write(iout,*) (non_zero_columns(j), j=1,n3d)
!  write(iout,*) (row_index(j), j=1,count)
!  write(iout,*) (packed_columns(j), j=1,count)
! 
! 
  WRITE(iout,2) n3d*(n3d+1)/2, number
1 FORMAT(/,10x,'Packed Matrix.  All Matrix Elements Held in Core')
2 FORMAT(/,10x,'Possible Non-Zero Elements =', i6,5x, 'Actual Non-Zero Elements =',i6)
END SUBROUTINE h_pack_lower_z
!***********************************************************************
!***********************************************************************
!*deck  @(#)pack_matrix_d
  subroutine pack_matrix_d(matrix,matrix_diagonal,packed_columns,non_zero_columns,row_index,n)
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
  INTEGER                           :: count
  count=0
  non_zero_columns(:) = 0
  write(iout,1)
  DO j=1,n
     DO i=j+1, n
        mat_el = matrix(i,j)
        IF(abs(mat_el) >= drop_tol) THEN 
           non_zero_columns(j) = non_zero_columns(j) + 1
           count = count + 1
           row_index(count) = i
           packed_columns(count) = mat_el
        END IF
     END DO
  END DO
  DO i=1,n
     matrix_diagonal(i) = matrix(i,i)
  END DO
  number = n*(n-1)/2
!  write(iout,*) non_zero_columns
!  write(iout,*) row_index(1:count)
!  write(iout,*) packed_columns(1:count)
!  write(iout,*) matrix_diagonal
  WRITE(iout,2) number, drop_tol, count
1 FORMAT(/,10x,'Packed Matrix.  All Matrix Elements Held in Core')
2 FORMAT(5x,'Total Number of Matrix Elements = ',i12,       &
       /,5x,'Effective Zero Cutoff           = ',d15.8,     &
       /,5x,'Number Exceeding Cutoff         = ',i12)
  END SUBROUTINE pack_matrix_d
!***********************************************************************
!***********************************************************************
!*deck  @(#)pack_matrix_z
  subroutine pack_matrix_z(matrix,matrix_diagonal,packed_columns,non_zero_columns,row_index,n)
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
  INTEGER                           :: count
  count=0
  non_zero_columns(:) = 0
  write(iout,1)
  DO j=1,n
     DO i=j+1, n
        mat_el = matrix(i,j)
        IF(abs(mat_el) >= drop_tol) THEN 
           non_zero_columns(j) = non_zero_columns(j) + 1
           count = count + 1
           row_index(count) = i
           packed_columns(count) = mat_el
        END IF
     END DO
  END DO
  DO i=1,n
     matrix_diagonal(i) = matrix(i,i)
  END DO
  number = n*(n+1)/2
  WRITE(iout,2) number, drop_tol, count
1 FORMAT(/,10x,'Packed Matrix.  All Matrix Elements Held in Core')
2 FORMAT(5x,'Total Number of Matrix Elements = ',i12,       &
       /,5x,'Effective Zero Cutoff           = ',d15.8,     &
       /,5x,'Number Exceeding Cutoff         = ',i12)
  END SUBROUTINE pack_matrix_z
!***********************************************************************
!***********************************************************************
  END  MODULE Pack_Matrix_Module
!***********************************************************************
!***********************************************************************
