!*********************************************************************
                      MODULE Coulomb_Variables_Module
                      USE input_output
!
                      IMPLICIT NONE
!
!                     CONSTANTS
!
  COMPLEX*16                                         :: eye               = (0.d0,1.d0)
  COMPLEX*16, EXTERNAL                               :: cgamma
  COMPLEX*16, EXTERNAL                               :: cpsi
  REAL*8                                             :: pi                = 3.1415926535897932384626433832795028841968D0
  REAL*8                                             :: sqrt2             = 1.4142135623730950488016887242096980785696d0
  REAL*8                                             :: invsqrt2          = .7071067811865475244008443621048490392848d0
  REAL*8                                             :: eulerc            = .577215664901532860606512090082402431042D0
  REAL*8                                             :: tm30              = 1.0D-30
  REAL*8                                             :: tm16              = 1.0D-16
  REAL*8                                             :: sqrt_2_div_pi     = .79788456080286535587989211986876373695173d0
  REAL*8                                             :: zero              =  0.d0
  REAL*8                                             :: quarter           = .25d0
  REAL*8                                             :: half              = .50d0
  REAL*8                                             :: one               = 1.0d0
  REAL*8                                             :: two               = 2.0d0
  REAL*8                                             :: three             = 3.0d0
  REAL*8                                             :: four              = 4.0d0
  REAL*8                                             :: five              = 5.0d0
  REAL*8                                             :: six               = 6.0d0
  REAL*8                                             :: seven             = 7.0d0
  REAL*8                                             :: eight             = 8.0d0
  REAL*8                                             :: nine              = 9.0d0
  REAL*8                                             :: ten               = 10.d0
  REAL*8                                             :: ten_2             = 1.0D2
  REAL*8                                             :: ten_4             = 1.0D4
  REAL*8                                             :: ten_6             = 1.0D6
  REAL*8                                             :: tol               = 1.D-14
  REAL                                               :: zero_r4           = 0.0E0
  REAL                                               :: half_r4           = 0.5E0
  REAL                                               :: one_r4            = 1.0E0
  REAL                                               :: six_r4            = 6.0E0
  REAL                                               :: ten_r4            = 1.0E1
  REAL                                               :: rl35              = 3.5E1
  REAL                                               :: aloge             = 0.4342945E0
  REAL*8                                             :: abort             = 2.0D+04
  REAL*8                                             :: abort2            = 4.0D4
  REAL*8, EXTERNAL                                   :: gamma
  REAL*8, EXTERNAL                                   :: c_arg
  INTEGER                                            :: int_zero          = 0
  INTEGER                                            :: int_one           = 1
  INTEGER                                            :: int_two           = 2
  INTEGER                                            :: int_three         = 3
  INTEGER                                            :: int_four          = 4
  INTEGER                                            :: int_five          = 5
  INTEGER                                            :: int_six           = 6
  INTEGER                                            :: int_seven         = 7
  INTEGER                                            :: int_eight         = 8
  INTEGER                                            :: int_nine          = 9
  INTEGER                                            :: int_ten           = 10
  INTEGER                                            :: iuo               = 70
!

!                     MAJOR VARIABLES
!
  REAL*8                                             :: energy
  REAL*8                                             :: charge
  REAL*8                                             :: k
  REAL*8                                             :: r_series 
  REAL*8                                             :: r_asymptotic
  LOGICAL                                            :: print_sigma_l
  LOGICAL                                            :: print_long_range_coefficients
  LOGICAL                                            :: print_convergence
  LOGICAL                                            :: print_short_range_coefficients
  INTEGER                                            :: series_size
  INTEGER                                            :: asymptotic_size
  REAL*8                                             :: r
  REAL*8                                             :: r_inv
  REAL*8                                             :: eta_in
  REAL*8                                             :: wronskian
  CHARACTER(LEN=3), EXTERNAL                         :: itoc
  CHARACTER(LEN=80)                                  :: title
  CHARACTER(LEN=8)                                   :: energy_class
  CHARACTER(LEN=1600)                                :: card
  CHARACTER(LEN=80)                                  :: cpass
  CHARACTER(LEN=80)                                  :: quantities_returned
  CHARACTER(LEN=80)                                  :: type
  INTEGER                                            :: angular_momentum
  INTEGER                                            :: l_max
  INTEGER                                            :: l_min
  INTEGER                                            :: number_of_r_values
  INTEGER                                            :: iexp
  INTEGER                                            :: ifail
  REAL*8                                             :: r_min
  REAL*8                                             :: r_max
  REAL*8                                             :: r_step
  REAL*8                                             :: fl
  REAL*8                                             :: dfl
  REAL*8                                             :: gl
  REAL*8                                             :: dgl
  REAL*8, DIMENSION(4)                               :: phase_factor
  REAL*8                                             :: u
  REAL*8                                             :: twou
  REAL*8                                             :: fouru

!
!                    ALLOCATABLES
!
  REAL*8, DIMENSION(:,:), ALLOCATABLE                :: fgd
  REAL*8, DIMENSION(:),   ALLOCATABLE                :: rho
  REAL*8, DIMENSION(:),   ALLOCATABLE                :: rho_inv
  TYPE series_coefficients
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: a
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: b
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: a_0
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: b_0
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: c_0
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: d_0
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: a_i
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: b_i
       REAL*8,   DIMENSION(:),                                       &
                 ALLOCATABLE                         :: e_0
  END TYPE series_coefficients
  TYPE(series_coefficients),                                         &
                  DIMENSION(:),                                      &
                  ALLOCATABLE                        :: power_series
  TYPE(series_coefficients),                                         &
                  DIMENSION(:),                                      &
                  ALLOCATABLE                        :: asymptotic_series
!====================================================================
!====================================================================
END MODULE Coulomb_Variables_Module
