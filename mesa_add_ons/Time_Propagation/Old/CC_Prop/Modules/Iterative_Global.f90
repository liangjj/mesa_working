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
  USE input_output
  USE prop_global
  IMPLICIT NONE
!
!                           Allocations for Real Variables
!
  REAL*8,       DIMENSION(:,:),     &
                ALLOCATABLE              :: vec_d
  REAL*8,       DIMENSION(:),       &
                ALLOCATABLE              :: h_vec_d
  REAL*8,       DIMENSION(:,:),     &
                ALLOCATABLE              :: s_mat_d
  REAL*8,       DIMENSION(:,:),     &
                ALLOCATABLE              :: s_mat_work_d
  REAL*8,       DIMENSION(:,:),     &
                ALLOCATABLE              :: h_mat_d
  REAL*8,       DIMENSION(:,:),     &
                ALLOCATABLE              :: h_mat_work_d 
  REAL*8,       DIMENSION(:),       &
                ALLOCATABLE              :: v_tri_d
  REAL*8,       DIMENSION(:),       &
                ALLOCATABLE              :: a, b
!
!                           Allocations for Complex Variables
  COMPLEX*16,   DIMENSION(:,:),     &
                ALLOCATABLE              :: vec_z
  COMPLEX*16,   DIMENSION(:),       &
                ALLOCATABLE              :: h_vec_z
  COMPLEX*16,   DIMENSION(:,:),     &
                ALLOCATABLE              :: s_mat_z
  COMPLEX*16,   DIMENSION(:,:),     &
                ALLOCATABLE              :: s_mat_work_z
  COMPLEX*16,   DIMENSION(:,:),     &
                ALLOCATABLE              :: h_mat_z
  COMPLEX*16,   DIMENSION(:,:),     &
                ALLOCATABLE              :: h_mat_work_z
  COMPLEX*16,   DIMENSION(:),       &
                ALLOCATABLE              :: v_tri_z
  REAL*8,       DIMENSION(:),       &
!
!
!                          Common variables to a number of routines
!
                ALLOCATABLE              :: eig
  COMPLEX*16,   DIMENSION(:,:),     &
                ALLOCATABLE              :: hamiltonian_z 
  COMPLEX*16,   DIMENSION(:,:),     &
                ALLOCATABLE              :: overlap_z 
  COMPLEX*16,   DIMENSION(:,:),     &
                ALLOCATABLE              :: eigenvectors_z
  REAL*8,       DIMENSION(:,:),     &
                ALLOCATABLE              :: hamiltonian_d 
  REAL*8,       DIMENSION(:,:),     &
                ALLOCATABLE              :: overlap_d 
  REAL*8,       DIMENSION(:,:),     &
                ALLOCATABLE              :: eigenvectors_d 
  REAL*8,       DIMENSION(:),       &
                ALLOCATABLE              :: eigenvalues
! 
!                          Allocations for scratch variables
!
  REAL*8,       DIMENSION(:),       &
                ALLOCATABLE              :: work_d
  COMPLEX*16,   DIMENSION(:),       &
                ALLOCATABLE              :: work_z
  REAL*8,       DIMENSION(:),       &
                ALLOCATABLE              :: rwork 
!
!                          Directives, Convergence and other variables
!
  LOGICAL                                :: non_orth
  REAL*8                                 :: error
  REAL*8                                 :: wfn_tst
  INTEGER                                :: lwork
  REAL*8                                 :: eps=1.d-08
  LOGICAL                                :: null_vec
  CHARACTER(LEN=16)                      :: orthogonalize, iterative_method
  CHARACTER(LEN=16)                      :: convergence
  CHARACTER(LEN=16)                      :: type
  INTEGER                                :: last_it, it_min
  INTEGER                                :: info
  REAL*8                                 :: anorm
!
!                          Keywords for Print Options
!
  CHARACTER (LEN=80), DIMENSION(10)      :: print_iterative
  LOGICAL,            DIMENSION(10)      :: log_iterative
  DATA print_iterative  / 'recursion_coefficients','iterative_vectors',       &
                        'h_on_vector','overlap_matrix', 'small_matrices',     &
                        'eigenvalues', 'eigenvectors','convergence_tests',    &
                        'solution', 'none' /
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Iterative_Global




