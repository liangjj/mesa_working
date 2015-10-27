MODULE Derived_Types
!***begin prologue     Derived_Types
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           time, finite element dvr, orthogonal polynomial
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            global shared variables for dvr library
!***description        this routine defines the global variables
!***                   and data needed for the dvrlib
!
!***references

!***routines called    
!***end prologue       Derived_Types
  USE accuracy
  USE Matrix_Print
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!---------------------------------------------------------------------
!                    Derived types and Allocated variables for 
!                             the dvr quantities
!---------------------------------------------------------------------
!
  TYPE(real_vector)                                  :: type_real_vector
  TYPE(real_matrix)                                  :: type_real_matrix
!
  TYPE even_matrices
    REAL(idp), DIMENSION(:,:),  ALLOCATABLE          :: tr
    REAL(idp), DIMENSION(:,:),  ALLOCATABLE          :: ham
  END TYPE even_matrices

  TYPE odd_matrices
    REAL(idp), DIMENSION(:,:),  ALLOCATABLE          :: tr
    REAL(idp), DIMENSION(:,:),  ALLOCATABLE          :: ham
  END TYPE odd_matrices

  TYPE final_matrices
    REAL(idp), DIMENSION(:,:),  ALLOCATABLE          :: tr
    REAL(idp), DIMENSION(:,:),  ALLOCATABLE          :: ham
  END TYPE final_matrices

  TYPE working_arrays
    TYPE (final_matrices),                                                   &
                     DIMENSION(:),    ALLOCATABLE    :: mat
    REAL(idp),       DIMENSION(:),    ALLOCATABLE    :: eigenvalues
    REAL(idp),       DIMENSION(:,:),  ALLOCATABLE    :: eigenvectors
    REAL(idp),       DIMENSION(:),    ALLOCATABLE    :: work
    REAL(idp),       DIMENSION(:),    ALLOCATABLE    :: lower_mat
    REAL(idp),       DIMENSION(:,:),  ALLOCATABLE    :: ham
    REAL(idp),       DIMENSION(:,:),  ALLOCATABLE    :: tr
  END TYPE working_arrays

  TYPE Regional
    REAL(idp),       DIMENSION(:),    ALLOCATABLE   :: q
    REAL(idp),       DIMENSION(:),    ALLOCATABLE   :: wt
    REAL(idp),       DIMENSION(:),    ALLOCATABLE   :: q_fac
    REAL(idp),       DIMENSION(:),    ALLOCATABLE   :: inv_q_fac
    REAL(idp),       DIMENSION(:),    ALLOCATABLE   :: inv_sqrt_q_fac
    REAL(idp),       DIMENSION(:,:),  ALLOCATABLE   :: p
    REAL(idp),       DIMENSION(:,:),  ALLOCATABLE   :: dp
    REAL(idp),       DIMENSION(:,:),  ALLOCATABLE   :: ddp
    REAL(idp),       DIMENSION(:,:),  ALLOCATABLE   :: bf
    REAL(idp),       DIMENSION(:,:),  ALLOCATABLE   :: dbf
    REAL(idp),       DIMENSION(:,:),  ALLOCATABLE   :: ddbf
    REAL(idp),       DIMENSION(:),    ALLOCATABLE   :: normalization
    REAL(idp),       DIMENSION(:),    ALLOCATABLE   :: edge
    REAL(idp),       DIMENSION(:),    ALLOCATABLE   :: pot
    TYPE (even_matrices)                            :: even
    TYPE (odd_matrices)                             :: odd
    TYPE (final_matrices),                                                    &
                     DIMENSION(:),    ALLOCATABLE   :: mat
    INTEGER                                         :: n_pts
    INTEGER                                         :: n_fun
    INTEGER                                         :: n_fixed
    INTEGER                                         :: first
    INTEGER                                         :: last
    INTEGER                                         :: fixed_point
    CHARACTER(LEN=8)                                :: type_quadrature
  END TYPE Regional

  TYPE cartesian
    CHARACTER(LEN=8)                                :: axis
  END TYPE cartesian

  TYPE spherical
    INTEGER                                         :: l_max
    INTEGER                                         :: m_max
    CHARACTER(LEN=8)                                :: axis
  END TYPE spherical

  TYPE spheroidal
    INTEGER                                         :: m_max
    REAL(idp)                                       :: R_ab
    REAL(idp)                                       :: Z_a
    REAL(idp)                                       :: Z_b
    CHARACTER(LEN=8)                                :: axis
  END TYPE spheroidal

  TYPE cylindrical
    INTEGER                                         :: m_max
    CHARACTER(LEN=8)                                :: axis
  END TYPE cylindrical

  TYPE potential
    REAL(idp)                                       :: depth
    INTEGER                                         :: nwell
    REAL(idp)                                       :: awell
    REAL(idp), DIMENSION(2)                         :: expnt
    REAL(idp), DIMENSION(2)                         :: amp
    INTEGER                                         :: n_p
    REAL(idp)                                       :: charge
    REAL(idp)                                       :: mass
    REAL(idp)                                       :: omega
    REAL(idp)                                       :: shift
    INTEGER                                         :: n_scale
    REAL(idp)                                       :: e_c
    REAL(idp)                                       :: dummy
    REAL(idp), DIMENSION(0:2)                       :: c
  END TYPE potential

  TYPE coordinates
    INTEGER                                         :: n_fix
    LOGICAL,           DIMENSION(2)                 :: drop
    LOGICAL,           DIMENSION(2)                 :: fix
    REAL(idp),         DIMENSION(:),   ALLOCATABLE  :: points
    REAL(idp),         DIMENSION(:),   ALLOCATABLE  :: weights
    TYPE(cartesian)                                 :: xyz
    TYPE(spherical)                                 :: r_theta
    TYPE(cylindrical)                               :: rho_z
    TYPE(spheroidal)                                :: xi_eta
    TYPE(lebedev)                                   :: theta_phi
    TYPE (regional),  DIMENSION(:),    ALLOCATABLE  :: reg
    TYPE (regional),  DIMENSION(:),    ALLOCATABLE  :: temp
  END TYPE coordinates

    TYPE (coordinates)                              :: grid
    TYPE (potential)                                :: pe
    TYPE (working_arrays)                           :: wa
!
!
END MODULE Derived_Types
