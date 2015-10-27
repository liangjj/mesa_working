!***********************************************************************
!**begin prologue     Atomic_Matrices
!**date written       080612   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time propagation
!**
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            This module packages all of the variables associated with
!***                  the use of atomic symmetry.
!***                  
!***description       
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!**references
!**modules needed     None directly
!**end prologue       Atomic_Matrices
!***********************************************************************
!***********************************************************************
                           MODULE Atomic_Matrices
!
                           IMPLICIT NONE
!
  TYPE state_labels
       INTEGER,   DIMENSION(:),             &
                  ALLOCATABLE                :: state_quantum_numbers
  END TYPE state_labels
!
  TYPE(state_labels),                       &
                  DIMENSION(:),             &
                  ALLOCATABLE                :: labels
  TYPE state_matrix
     REAL*8,      DIMENSION(:),             &
                  ALLOCATABLE               :: state_h_matrix
     REAL*8,      DIMENSION(:),             &
                  ALLOCATABLE               :: state_s_matrix
     REAL*8,      DIMENSION(:),             &
                  ALLOCATABLE               :: state_eigen_values
     REAL*8,      DIMENSION(:,:),           &
                  ALLOCATABLE               :: state_eigen_vectors
     INTEGER,     DIMENSION(:),             &
                  ALLOCATABLE               :: h_non_zero_columns
     INTEGER,     DIMENSION(:),             &
                  ALLOCATABLE               :: h_row_index
     REAL*8,      DIMENSION(:),             &
                  ALLOCATABLE               :: h_packed_columns
     REAL*8,      DIMENSION(:),             &
                  ALLOCATABLE               :: h_diagonal
     INTEGER                                :: h_number
     INTEGER,     DIMENSION(:),             &
                  ALLOCATABLE               :: s_non_zero_columns
     INTEGER,     DIMENSION(:),             &
                  ALLOCATABLE               :: s_row_index
     REAL*8,      DIMENSION(:),             &
                  ALLOCATABLE               :: s_packed_columns
     REAL*8,      DIMENSION(:),             &
                  ALLOCATABLE               :: s_diagonal
     INTEGER                                :: s_number
     INTEGER,     DIMENSION(:),             &
                  ALLOCATABLE               :: eig_non_zero_columns
     INTEGER,     DIMENSION(:),             &
                  ALLOCATABLE               :: eig_row_index
     REAL*8,      DIMENSION(:),             &
                  ALLOCATABLE               :: eig_packed_columns
     REAL*8,      DIMENSION(:),             &
                  ALLOCATABLE               :: eig_diagonal
     INTEGER                                :: eig_number
     INTEGER,   DIMENSION(:),               &
                  ALLOCATABLE               :: matrix_pointer
  END TYPE state_matrix
!
  TYPE (state_matrix),                      &
                  DIMENSION(:),             &
                  ALLOCATABLE               :: state_mat
!
  TYPE dipole_matrix
     REAL*8,      DIMENSION(:,:),           &
                  ALLOCATABLE               :: dipole_coupling_matrix
     INTEGER,     DIMENSION(:),             &
                  ALLOCATABLE               :: dipole_non_zero_columns
     INTEGER,     DIMENSION(:),             &
                  ALLOCATABLE               :: dipole_row_index
     REAL*8,      DIMENSION(:),             &
                  ALLOCATABLE               :: dipole_packed_columns
     INTEGER                                :: dipole_number
  END TYPE dipole_matrix
!
  TYPE (dipole_matrix),                     &
                  DIMENSION(:),             &
                  ALLOCATABLE               :: dipole_mat
!
     INTEGER,   DIMENSION(:),               &
                  ALLOCATABLE               :: matrix_size
!
    INTEGER,   DIMENSION(:),                &
                  ALLOCATABLE               :: state_matrix_size
    INTEGER,   DIMENSION(:),                &
                  ALLOCATABLE               :: final_state_matrix_size
!
  TYPE channel_labels
       INTEGER,   DIMENSION(:),             &
                  ALLOCATABLE                :: channel_quantum_numbers
  END TYPE channel_labels
!
  TYPE(channel_labels),                     &
                  DIMENSION(:),             &
                  ALLOCATABLE                :: channel
!
  TYPE channel_matrix
     REAL*8,      DIMENSION(:,:),           &
                  ALLOCATABLE               :: channel_h_matrix_d
     REAL*8,      DIMENSION(:,:),           &
                  ALLOCATABLE               :: channel_s_matrix_d
     COMPLEX*16,  DIMENSION(:,:),           &
                  ALLOCATABLE               :: channel_h_matrix_z
     COMPLEX*16,  DIMENSION(:,:),           &
                  ALLOCATABLE               :: channel_s_matrix_z
  END TYPE channel_matrix
!
  TYPE (channel_matrix),                    &
                  DIMENSION(:,:),           &
                  ALLOCATABLE               :: channel_mat
!
  INTEGER,        DIMENSION(:,:),           &
                  ALLOCATABLE               :: channel_buf
!
  REAL*8,         DIMENSION(:),             &
                  ALLOCATABLE               :: atomic_eigen_values
!
  INTEGER                                   :: L_Max
  INTEGER                                   :: locate_L
  INTEGER                                   :: locate_L_plus_one
  INTEGER                                   :: locate_L_minus_one
  INTEGER                                   :: number_of_channels
  INTEGER                                   :: state_tri_size
  INTEGER                                   :: channel_size
  INTEGER                                   :: final_total_matrix_size
!***********************************************************************
!***********************************************************************
  END  MODULE Atomic_Matrices
!***********************************************************************
!***********************************************************************
