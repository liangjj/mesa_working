!***********************************************************************
!**begin prologue     Global_Time_Propagation_Module
!**date written       080612   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time propagation
!**
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            This module packages all of the variables needed directly
!***                  in the propagation code. This includes arrays as well as scalars.
!***                  It does not include variables which are associated with the
!***                  iterative and packing modules.
!***                  
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
!**end prologue       Time_Propagation_Input_Data_Module
!***********************************************************************
!***********************************************************************
                           MODULE Global_Time_Propagation_Module
!
                           IMPLICIT NONE
!
!
  REAL*8,         DIMENSION(:),             &
                  ALLOCATABLE               :: s_vec_in_d
  REAL*8,         DIMENSION(:),             &
                  ALLOCATABLE               :: s_vec_out_d
  COMPLEX*16,     DIMENSION(:),             &
                  ALLOCATABLE               :: s_vec_in_z
  COMPLEX*16,     DIMENSION(:),             &
                  ALLOCATABLE               :: s_vec_out_z
!
!
     REAL*8,      DIMENSION(:,:),           &
                  ALLOCATABLE               :: scratch_matrix
!
     REAL*8,      DIMENSION(:),             &
                  ALLOCATABLE               :: scratch_tri
!
  INTEGER                                   :: L
  INTEGER                                   :: L_Prime
  INTEGER                                   :: number_of_splines
  INTEGER                                   :: number_of_correlation_terms
  INTEGER                                   :: spline_order
  INTEGER                                   :: number_of_splines_i
  INTEGER                                   :: number_of_splines_j
  INTEGER                                   :: first_spline
  INTEGER                                   :: last_spline
  INTEGER, DIMENSION(:), ALLOCATABLE        :: spline_array_i
  INTEGER, DIMENSION(:), ALLOCATABLE        :: spline_array_j
  INTEGER                                   :: ic
  INTEGER                                   :: jc
  INTEGER                                   :: is
  INTEGER                                   :: js
  INTEGER                                   :: ii
  INTEGER                                   :: jj
  INTEGER                                   :: new_ind
  INTEGER                                   :: old_ind
  INTEGER                                   :: old_ij
  INTEGER                                   :: new_ij
  INTEGER                                   :: first
  INTEGER                                   :: last
  INTEGER                                   :: first_matrix_index
  INTEGER                                   :: last_matrix_index
  INTEGER                                   :: first_new_matrix_index
  INTEGER                                   :: last_new_matrix_index
  INTEGER                                   :: tri_size
  INTEGER                                   :: new_size
  INTEGER                                   :: n_dim
  INTEGER                                   :: new_tri_size
  INTEGER                                   :: count
  INTEGER                                   :: final_total_matrix_size
  INTEGER                                   :: old_matrix_count_i
  INTEGER                                   :: old_matrix_count_j
  INTEGER                                   :: new_matrix_count_i
  INTEGER                                   :: new_matrix_count_j
  INTEGER                                   :: non_zero_overlap_elements
  INTEGER                                   :: non_zero_cholesky_elements
  INTEGER                                   :: non_zero_hamiltonian_elements
  INTEGER                                   :: lenbuf
  INTEGER                                   :: remove
  LOGICAL                                   :: channel_format
  LOGICAL, DIMENSION(2)                     :: spline_removal
  LOGICAL                                   :: remove_last_spline
  LOGICAL                                   :: remove_next_to_last_spline
  LOGICAL                                   :: print_channel_matrices
  LOGICAL                                   :: print_cc
  LOGICAL                                   :: print_packed_matrices
  LOGICAL                                   :: print_input_matrices
  LOGICAL                                   :: print_internal_matrices
  LOGICAL                                   :: print_buffers
  LOGICAL                                   :: print_packed_disk_buffers
  LOGICAL                                   :: channel_labels_only
  LOGICAL                                   :: dipole_matrices
  LOGICAL                                   :: write_pointers
  LOGICAL                                   :: write_channel_labels
  LOGICAL                                   :: ifeig=.false.
  LOGICAL                                   :: diagonalize_only
  LOGICAL                                   :: diagonalize_matrix
  LOGICAL                                   :: read_only
  LOGICAL                                   :: reformat_input_matrix_only
  LOGICAL                                   :: test_s_inverse
  LOGICAL                                   :: cholesky_from_disk
  LOGICAL                                   :: non_orth
  LOGICAL                                   :: packed_matrices
  LOGICAL                                   :: in_core
  LOGICAL                                   :: full_cholesky_to_disk
  LOGICAL                                   :: overlap_to_disk
  LOGICAL                                   :: time_dependent_potential
  LOGICAL                                   :: to_standard_eigenvalue_problem
  CHARACTER (LEN=1600)                      :: card
  CHARACTER (LEN=80)                        :: cpass
  CHARACTER (LEN=80)                        :: title
  CHARACTER (LEN=80)                        :: ham_file
  CHARACTER (LEN=80)                        :: file_format
  CHARACTER (LEN=80)                        :: input_matrices
  CHARACTER (LEN=80)                        :: output_matrices
  CHARACTER (LEN=128)                       :: cholesky_file_name
  CHARACTER (LEN=128)                       :: packed_file_name
  CHARACTER (LEN=128)                       :: state_file_name
  CHARACTER (LEN=80)                        :: reformatting_control
  CHARACTER (LEN=80)                        :: iteration_method
  CHARACTER (LEN=8)                         :: input_output_vector
  CHARACTER (LEN=80)                        :: matrix_type
  CHARACTER (LEN=80)                        :: matrix_file
  CHARACTER (LEN=80)                        :: matrix_directive
  CHARACTER (LEN=3)                         :: L_1
  CHARACTER (LEN=3)                         :: L_2
  CHARACTER (LEN=4)                         :: parity
  CHARACTER (LEN=16)                        :: type_pulse
  CHARACTER (LEN=8)                         :: use_atomic_symmetry
  CHARACTER (LEN=8)                         :: file_key
  CHARACTER (LEN=16)                        :: peak_type = 'intensity'
  CHARACTER (LEN=16)                        :: spatial_representation = 'spline'
  CHARACTER (LEN=16)                        :: species
  REAL*8                                    :: smallest
  REAL*8                                    :: largest
  REAL*8                                    :: drop_overlap
  REAL*8                                    :: drop_hamiltonian
  REAL*8                                    :: drop_cholesky
  REAL*8                                    :: drop_dipole
  REAL*8                                    :: drop_vectors
  REAL*8                                    :: electric_field
  REAL*8                                    :: peak_electric_field=.25d-04
  REAL*8                                    :: peak_intensity
  REAL*8                                    :: pulse_duration = 30.d0
  REAL*8                                    :: ramp_time = 5.d0
  REAL*8                                    :: photon_energy
  REAL*8                                    :: time_one_optical_cycle
  REAL*8                                    :: energy_cutoff = 10.d0
  INTEGER                                   :: number_optical_cycles
  INTEGER                                   :: steps_per_optical_cycle
!***********************************************************************
!***********************************************************************
  END  MODULE Global_Time_Propagation_Module
!***********************************************************************
!***********************************************************************
