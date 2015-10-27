!
MODULE FEDVR_Derived_Types
!***begin prologue     FEDVR_Derived_Types
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
!***end prologue       FEDVR_Derived_Types
  USE accuracy
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
  TYPE dvr_matrices
    REAL(idp), DIMENSION(:,:), ALLOCATABLE       :: ham
    REAL(idp), DIMENSION(:,:), ALLOCATABLE       :: eigenvectors
    REAL(idp), DIMENSION(:),   ALLOCATABLE       :: eigenvalues
    REAL(idp), DIMENSION(:,:), ALLOCATABLE       :: rhs
    REAL(idp), DIMENSION(:,:), ALLOCATABLE       :: computed_solution
    REAL(idp), DIMENSION(:,:), ALLOCATABLE       :: exact_solution
    REAL(idp), DIMENSION(:),   ALLOCATABLE       :: work
    REAL(idp), DIMENSION(:),   ALLOCATABLE       :: points
    REAL(idp), DIMENSION(:),   ALLOCATABLE       :: lower
    INTEGER,   DIMENSION(:),   ALLOCATABLE       :: ipvt
    INTEGER                                      :: number_of_right_hand_sides
    CHARACTER(LEN=80)                            :: type_inhomo
    LOGICAL                                      :: poisson
  END TYPE dvr_matrices

  TYPE pt_wt
       REAL(idp), DIMENSION(:), ALLOCATABLE      :: qr
       REAL(idp), DIMENSION(:), ALLOCATABLE      :: wtr
       REAL(idp), DIMENSION(:), ALLOCATABLE      :: qr_fac
       REAL(idp), DIMENSION(:), ALLOCATABLE      :: inv_sqrt_qr_fac
       REAL(idp), DIMENSION(:), ALLOCATABLE      :: inv_qr_fac
  END TYPE pt_wt

  TYPE functions
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: pr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: dpr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ddpr
       REAL(idp), DIMENSION(:),   ALLOCATABLE    :: normalization
  END TYPE functions

  TYPE odd_functions
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: pr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: dpr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ddpr
  END TYPE odd_functions

  TYPE fourier_functions
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: pr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: dpr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ddpr
       REAL(idp), DIMENSION(:),   ALLOCATABLE    :: normalization
  END TYPE fourier_functions

  TYPE matrices
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE matrices

  TYPE fourier_matrices
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE fourier_matrices

  TYPE laguerre_matrices
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE laguerre_matrices

  TYPE hermite_matrices
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE hermite_matrices

  TYPE odd_matrices
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE odd_matrices

  TYPE hamiltonian
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: tr
       REAL(idp), DIMENSION(:,:), ALLOCATABLE    :: ham
  END TYPE hamiltonian

  TYPE vectors
       REAL(idp), DIMENSION(:),   ALLOCATABLE    :: vr
       REAL(idp), DIMENSION(:),   ALLOCATABLE    :: s
  END TYPE vectors

  TYPE well
       REAL(idp)                                 :: depth
       INTEGER                                   :: nwell
       REAL(idp)                                 :: awell
  END TYPE well

  TYPE exponential
       REAL(idp), DIMENSION(2)                   :: expnt
       REAL(idp), DIMENSION(2)                   :: amp
  END TYPE exponential

  TYPE sum_exponential
       REAL(idp), DIMENSION(2)                   :: expnt
       REAL(idp), DIMENSION(2)                   :: amp
  END TYPE sum_exponential

  TYPE power_exponential
       REAL(idp), DIMENSION(2)                   :: expnt
       REAL(idp), DIMENSION(2)                   :: amp
       INTEGER                                   :: n_p
  END TYPE power_exponential

  TYPE yukawa
       REAL(idp), DIMENSION(2)                   :: expnt
       REAL(idp), DIMENSION(2)                   :: amp
  END TYPE yukawa

  TYPE coulomb
       REAL(idp)                                 :: charge
  END TYPE coulomb

  TYPE eberlonium
       REAL(idp)                                 :: charge
       REAL(idp), DIMENSION(2)                   :: amp
       INTEGER                                   :: n_p
  END TYPE eberlonium

  TYPE harmonic_oscillator
       REAL(idp)                                 :: mass
       REAL(idp)                                 :: omega
  END TYPE harmonic_oscillator

  TYPE expres
       REAL(idp), DIMENSION(2)                   :: expnt
       REAL(idp), DIMENSION(2)                   :: amp
       REAL(idp)                                 :: shift
  END TYPE expres

  TYPE periodic
       INTEGER                                   :: n_scale
       REAL(idp)                                 :: e_c
       REAL(idp)                                 :: awell
  END TYPE periodic

  TYPE inverse_r4
       REAL(idp)                                 :: dummy
  END TYPE inverse_r4 

  TYPE potential
       TYPE (well)                               :: well
       TYPE (yukawa)                             :: yukawa
       TYPE (exponential)                        :: exponential
       TYPE (coulomb)                            :: coulomb
       TYPE (expres)                             :: expres
       TYPE (harmonic_oscillator)                :: harmonic_oscillator
       TYPE (eberlonium)                         :: eberlonium
       TYPE (power_exponential)                  :: power_exp
       TYPE (sum_exponential)                    :: sum_exp
       TYPE (inverse_r4)                         :: inv_r4
       TYPE (periodic)                           :: periodic
       CHARACTER(LEN=80)                         :: type
       REAL(idp), DIMENSION(:),   ALLOCATABLE    :: vec
  END TYPE potential

  TYPE coordinates
       CHARACTER(LEN=24)                         :: label
       INTEGER                                   :: num_reg
       INTEGER, DIMENSION(:),     ALLOCATABLE    :: num_pts_reg
       LOGICAL, DIMENSION(2)                     :: drop_pt
       LOGICAL, DIMENSION(2)                     :: fix_pt
       INTEGER                                   :: num_fixed
       REAL(idp), DIMENSION(:),   ALLOCATABLE    :: grid_points
       REAL(idp), DIMENSION(:),   ALLOCATABLE    :: grid_weights
       TYPE(pt_wt),               DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_pt_wt
       TYPE(functions),           DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_poly
       TYPE(odd_functions),       DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_poly_odd
       TYPE(fourier_functions),   DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_poly_fourier
       TYPE(matrices),            DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_mat
       TYPE(odd_matrices),        DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_mat_odd
       TYPE(fourier_matrices),    DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_mat_fourier
       TYPE(hermite_matrices),    DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_mat_hermite
       TYPE(laguerre_matrices),   DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_mat_laguerre
       TYPE(vectors),             DIMENSION(:),                                 &
                                  ALLOCATABLE    :: reg_vec
       TYPE(hamiltonian),         DIMENSION(:,:),                               &
                                  ALLOCATABLE    :: reg_ham
  END TYPE coordinates

  TYPE (coordinates),           DIMENSION(:),                                   &
                                ALLOCATABLE      :: reg_grid

  TYPE (potential),             DIMENSION(:),                                   &
                                ALLOCATABLE      :: reg_pot

  TYPE(dvr_matrices),           DIMENSION(:),                                   &
                                ALLOCATABLE      :: dvr_mat

  TYPE buffer
    REAL(idp),    DIMENSION(:), ALLOCATABLE                                     &
                                                 :: d, hbuf
    REAL(idp),    DIMENSION(:,:), ALLOCATABLE                                   &
                                                 :: hibuf
  END TYPE buffer

  TYPE(buffer), DIMENSION(:),                                                   &
                ALLOCATABLE                      :: buf    
  REAL(idp), DIMENSION(:), ALLOCATABLE           :: v_pert
!
!
END MODULE FEDVR_Derived_Types
