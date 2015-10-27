  MODULE lanczos_global
!deck lanczos_global.f
!**begin prologue     lanczos_global
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            global arrays for time propagator.
!**references
!**routines called
!**end prologue       lanczos_global
  USE input_output
  USE prop_global
  IMPLICIT NONE
!
!
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: vec_d
  REAL*8,       DIMENSION(:),                 &          
                ALLOCATABLE                     :: h_vec_d
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: s_vec_d
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: s_mat_d
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: s_mat_work_d
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: h_small_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: v_tri_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: work_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: vscr_d
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: rwork 
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: eig
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: a, b
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: hamiltonian_d 
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: overlap_d 
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: eigenvectors_d 
  REAL*8,       DIMENSION(:,:),               &
                ALLOCATABLE                     :: eigen_vectors 
  REAL*8,       DIMENSION(:),                 &
                ALLOCATABLE                     :: eigenvalues
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: vec_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: h_vec_z
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: s_vec_z
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: s_mat_z
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: s_mat_work_z
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: h_small_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: v_tri_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: work_z
  COMPLEX*16,   DIMENSION(:),                 &
                ALLOCATABLE                     :: vscr_z
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: hamiltonian_z 
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: overlap_z 
  COMPLEX*16,   DIMENSION(:,:),               &
                ALLOCATABLE                     :: eigenvectors_z
! 
  LOGICAL                                       :: non_orth
  REAL*8                                        :: error
  REAL*8                                        :: wfn_tst
  INTEGER                                       :: lwork
  INTEGER                                       :: maximum_number_of_time_subintervals
  REAL*8                                        :: eps=1.d-08
  LOGICAL                                       :: null_vec
  LOGICAL                                       :: schmidt
  CHARACTER(LEN=16)                             :: convergence
  CHARACTER(LEN=16)                             :: lanczos_convergence
  CHARACTER(LEN=16)                             :: type
  INTEGER                                       :: last_it
  INTEGER                                       :: info
  REAL*8                                        :: anorm
  REAL*8                                        :: norm
  INTEGER                                       :: total_number_of_iterations
  REAL*8                                        :: t_start
  REAL*8                                        :: t_end
  REAL*8                                        :: save_deltat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE lanczos_global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




