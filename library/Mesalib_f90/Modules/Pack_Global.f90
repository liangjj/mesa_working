!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE Pack_Global
!deck Pack_Global.f
!**begin prologue     Pack_Global
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            global arrays for packing routines.
!**references
!**routines called
!**end prologue       Pack_Global
  IMPLICIT NONE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE matrix
!                           Allocations for Real Variables
!
     INTEGER,      DIMENSION(:),            &
                   ALLOCATABLE                :: col_index  
     INTEGER,      DIMENSION(:),            &
                   ALLOCATABLE                :: row_index  
     INTEGER,      DIMENSION(:),            &
                   ALLOCATABLE                :: col_index_v  
     INTEGER,      DIMENSION(:),            &
                   ALLOCATABLE                :: row_index_v  
     INTEGER                                &
                                              :: number
     INTEGER,      DIMENSION(:),            &
                   ALLOCATABLE                :: non_zero_rows
     INTEGER,      DIMENSION(:),            &
                   ALLOCATABLE                :: non_zero_columns
     INTEGER,      DIMENSION(:),            &
                   ALLOCATABLE                :: non_zero_rows_v
     INTEGER,      DIMENSION(:),            &
                   ALLOCATABLE                :: non_zero_columns_v
     INTEGER                                &
                                              :: max_col
     INTEGER                                &
                                              :: UNIT_NUMBER
     CHARACTER (LEN=16)                     &
                                              :: UNIT_NAME
     REAL*8,       DIMENSION(:),            &
                   ALLOCATABLE                :: packed_matrix_d
     REAL*8,       DIMENSION(:),            &
                   ALLOCATABLE                :: packed_rows_d
     REAL*8,       DIMENSION(:),            &
                   ALLOCATABLE                :: packed_columns_d
     REAL*8,       DIMENSION(:),            &
                   ALLOCATABLE                :: packed_rows_v_d
     REAL*8,       DIMENSION(:),            &
                   ALLOCATABLE                :: packed_columns_v_d
     REAL*8,       DIMENSION(:),            &
                   ALLOCATABLE                :: matrix_diagonal_d
!
!                           Allocations for Complex Variables
     COMPLEX *16,  DIMENSION(:),            &
                   ALLOCATABLE                :: matrix_diagonal_z
     COMPLEX*16,   DIMENSION(:),            &
                   ALLOCATABLE                :: packed_matrix_z
     COMPLEX*16,   DIMENSION(:),            &
                   ALLOCATABLE                :: packed_rows_z
     COMPLEX*16,   DIMENSION(:),            &
                   ALLOCATABLE                :: packed_columns_z
     COMPLEX*16,   DIMENSION(:),            &
                   ALLOCATABLE                :: packed_rows_v_z
     COMPLEX*16,   DIMENSION(:),            &
                   ALLOCATABLE                :: packed_columns_v_z
!
  END TYPE matrix
  TYPE (matrix),   DIMENSION(:),            &
                   ALLOCATABLE                :: mat_var
  TYPE (matrix),   DIMENSION(:,:),          &
                   ALLOCATABLE                :: chan_mat
! 
!
!                          Directives, Convergence and other variables
!
  LOGICAL                                     :: print_packed
  INTEGER                                     :: pointer
  INTEGER                                     :: trips
  INTEGER                                     :: left_over
  LOGICAL                                     :: packed
  REAL*8                                      :: drop_tol
  CHARACTER(LEN=24)                           :: type_packing
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Pack_Global




