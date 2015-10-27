!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE Iterative_Global
!deck Iterative_Global.f
!**begin prologue     Iterative_Global
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            global arrays for time propagator.
!**references
!**routines called
!**end prologue       Iterative_Global
  IMPLICIT NONE
!
  INTEGER,      DIMENSION(:,:),     &
                ALLOCATABLE              :: ibuf
!                           Allocations for Real Variables
!
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: v_in_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: v_out_d
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: vec_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: h_vec_d
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: s_vec_d
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: h_vectors_d
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: s_mat_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: s_mat_tri_d
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: s_mat_work_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: s_mat_work_tri_d
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: h_mat_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: sub_diagonal
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: h_mat_tri_d
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: h_mat_work_d 
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: h_mat_work_tri_d 
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: h_buf_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: v_tri_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: lanczos_tri_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: diag_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: rhs_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: small_rhs_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: small_rhs_work_d
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: residual_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: rhs_tri_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: soln_tri_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: a, b
!
!                           Allocations for Complex Variables
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: v_in_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: v_out_z
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: vec_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: h_vec_z
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: s_vec_z
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: h_vectors_z
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: s_mat_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: s_mat_tri_z
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: s_mat_work_z
  COMPLEX*16,   DIMENSION(:)  ,               &
                ALLOCATABLE                     :: s_mat_work_tri_z
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: h_mat_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: h_mat_tri_z
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: h_mat_work_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: h_mat_work_tri_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: h_buf_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: v_tri_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: lanczos_tri_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: diag_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: rhs_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: small_rhs_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: small_rhs_work_z
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: residual_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: rhs_tri_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: soln_tri_z
!
!
!                          Common variables to a number of routines
!
!                           Allocations for Real Variables
!
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: eig
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: hamiltonian_d 
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: matrix_d 
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: overlap_d 
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: triangle_hamiltonian_d 
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: triangle_overlap_d 
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: triangle_matrix_d 
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: upper_d 
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: lower_d 
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: eigenvectors_d 
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: eigen_vectors 
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: eigen_values
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: working_eigen_values
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: eig_previous
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: eigen_vectors_previous 
  REAL,         DIMENSION(:),                 &
                ALLOCATABLE                     :: total_time
  REAL                                        &           :: local_time
!
!                           Allocations for complex Variables
!
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: hamiltonian_z 
 COMPLEX*16,   DIMENSION(:,:),                &
                ALLOCATABLE                     :: matrix_z 
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: overlap_z 
  COMPLEX*16,   DIMENSION(:)  ,               &
                ALLOCATABLE                     :: triangle_hamiltonian_z 
  COMPLEX*16,   DIMENSION(:)  ,               &
                ALLOCATABLE                     :: triangle_matrix_z 
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: triangle_overlap_z 
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: upper_z 
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: lower_z 
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: eigenvectors_z
! 
!                          Allocations for scratch variables
!
  INTEGER,      DIMENSION(:),                 &
                ALLOCATABLE                     :: ipvt
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: work_d
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: work_z
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: rwork 
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: local_scratch_d
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: local_scratch_z
!
!                          Directives, Convergence and other variables
!
  REAL*8                                        :: error
  REAL*8                                        :: wfn_tst
  REAL*8                                        :: eps=1.d-10
  REAL*8                                        :: near_zero = 1.d-08
  INTEGER                                       :: lwork
  LOGICAL                                       :: null_vec
  LOGICAL                                       :: triangle
  LOGICAL                                       :: compute_energy
  CHARACTER(LEN=16)                             :: orthogonalize
  CHARACTER(LEN=24)                             :: iterative_method
  CHARACTER(LEN=16)                             :: convergence
  CHARACTER(LEN=16)                             :: type
  INTEGER                                       :: last_it
  INTEGER                                       :: it_min
  INTEGER                                       :: non_zero
  INTEGER                                       :: info
  REAL*8                                        :: anorm
  REAL*8                                        :: norm
!
!                          Keywords for Print Options
!
  CHARACTER (LEN=80), DIMENSION(10)             :: print_iterative
  LOGICAL,            DIMENSION(10)             :: log_iterative
  DATA print_iterative  / 'recursion_coefficients','iterative_vectors',                 &
                        'h_on_vector','overlap_matrix', 'small_matrices',               &
                        'eigenvalues', 'eigenvectors','convergence_tests',              &
                        'solution', 'none' /
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Iterative_Global




