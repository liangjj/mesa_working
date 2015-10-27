!
MODULE Grid_Defined_Types
!***begin prologue     Grid_Defined_Types
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source             
!***purpose            grid defined type
!***description        
!***                   
!
!***references

!***routines called    
!***end prologue       Grid_Defined_Types
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
!                             grid program
!---------------------------------------------------------------------
!
!
!                There is a type call atoms.  That type has all the information on the
!                geometry and charges of the atomic centers, including any dummy center associated
!                with the scattering.
!                There is a sub-type of atoms, called shellsl.  Shells has information on how the atomic
!                space is divided into shells.  Each shell has radial and angular coordinates
!                described in the type coodinates, The coordinates have a sub-type regions.
!                This contains all the information for that coordinate in a given region.  Thus
!                all the important information is inherited as atoms%shl%coord%reg.
!
  INTEGER                                              :: n_fixed
  LOGICAL,            DIMENSION(2)                     :: fixed
  LOGICAL,            DIMENSION(2)                     :: fix
  LOGICAL,            DIMENSION(2)                     :: drop
  REAL(idp),          DIMENSION(2)                     :: edge
  INTEGER                                              :: nply
  INTEGER                                              :: nsubr
  INTEGER                                              :: nblock
  INTEGER,             PARAMETER                       :: maxreg=5000
  LOGICAL                                              :: automte
  LOGICAL                                              :: defaults
  LOGICAL                                              :: reuse_sector_information
  REAL(idp)                                            :: step
  REAL(idp)                                            :: boundl
  REAL(idp)                                            :: boundr
  REAL(isp)                                            :: t_start
  REAL(isp)                                            :: t_end
  TYPE(real_vector)                                    :: type_real_vector
  TYPE(real_matrix)                                    :: type_real_matrix
  REAL(idp),            DIMENSION(:,:), ALLOCATABLE    :: full_cartesian_grid
  REAL(idp),            DIMENSION(:,:), ALLOCATABLE    :: full_spherical_grid
  REAL(idp),            DIMENSION(:)  , ALLOCATABLE    :: full_weights
  REAL(idp),            DIMENSION(:)  , ALLOCATABLE    :: scr
!
  TYPE FINAL_MATRICES
    REAL(idp), DIMENSION(:,:),  ALLOCATABLE            :: tr
    REAL(idp), DIMENSION(:,:),  ALLOCATABLE            :: ham
  END TYPE FINAL_MATRICES
!
  TYPE EVEN_MATRICES
    REAL(idp), DIMENSION(:,:),  ALLOCATABLE            :: tr
    REAL(idp), DIMENSION(:,:),  ALLOCATABLE            :: ham
  END TYPE EVEN_MATRICES

  TYPE ODD_MATRICES
    REAL(idp), DIMENSION(:,:),  ALLOCATABLE            :: tr
    REAL(idp), DIMENSION(:,:),  ALLOCATABLE            :: ham
  END TYPE ODD_MATRICES
!
  TYPE REGIONAL
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: q
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: q_fac
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: inv_q_fac
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: inv_sqrt_q_fac
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: p
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: dp
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: ddp
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: bf
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: dbf
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: ddbf
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: normalization
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: wt
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: theta
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: phi
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: cos_theta
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: sin_theta
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: cos_phi
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: sin_phi
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: pot
    REAL(idp),          DIMENSION(2)                   :: edge
    REAL(idp)                                          :: start
    REAL(idp)                                          :: end
    CHARACTER(LEN=16)                                  :: type_quadrature
    REAL(idp)                                          :: increment
    INTEGER                                            :: n_pts
    INTEGER                                            :: n_fixed
    INTEGER                                            :: first
    INTEGER                                            :: last
    INTEGER                                            :: n_fun
    INTEGER                                            :: fixed_point
    TYPE (FINAL_MATRICES), DIMENSION(:), ALLOCATABLE   :: mat
    TYPE (EVEN_MATRICES)                               :: even
    TYPE (ODD_MATRICES)                                :: odd
  END TYPE REGIONAL
!
  TYPE BOUNDARY_CONDITIONS
    INTEGER                                            :: n_fixed
    LOGICAL,            DIMENSION(2)                   :: fixed
    LOGICAL,            DIMENSION(2)                   :: fix
    LOGICAL,            DIMENSION(2)                   :: drop
    REAL(idp),          DIMENSION(2)                   :: edge
  END TYPE BOUNDARY_CONDITIONS
!
  TYPE LEBEDEV
    REAL(idp),          DIMENSION(:,:),  ALLOCATABLE   :: q
    REAL(idp),          DIMENSION(:,:),  ALLOCATABLE   :: pts
    REAL(idp)                                          :: wt
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: theta
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: phi
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: cos_theta
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: sin_theta
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: cos_phi
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: sin_phi
    REAL(idp)                                          :: a
    REAL(idp)                                          :: b
    REAL(idp)                                          :: c
    REAL(idp)                                          :: sumwt
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: w
    INTEGER                                            :: n_pts
    INTEGER                                            :: nphi
    INTEGER                                            :: lmax
    INTEGER                                            :: mmax
    INTEGER                                            :: code
    INTEGER                                            :: num
    INTEGER                                            :: nleb
    INTEGER                                            :: nang
    INTEGER                                            :: n_phy
  END TYPE LEBEDEV
!
  TYPE THETA
    TYPE (REGIONAL),    DIMENSION(:),    ALLOCATABLE   :: reg
    TYPE(BOUNDARY_CONDITIONS)                          :: bc
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: q
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: q_fac
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: inv_q_fac
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: inv_sqrt_q_fac
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: p
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: dp
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: ddp
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: bf
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: dbf
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: ddbf
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: normalization
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: wt
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: theta
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: cos_theta
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: sin_theta
    REAL(idp)                                          :: start
    REAL(idp)                                          :: end
    REAL(idp),          DIMENSION(2)                   :: edge
    CHARACTER(LEN=16)                                  :: type_quadrature
    REAL(idp)                                          :: increment
    INTEGER                                            :: fixed_point
    INTEGER                                            :: n_phy
    INTEGER                                            :: n_pts
    INTEGER                                            :: n_fun
    INTEGER                                            :: n_reg
  END TYPE THETA
!
  TYPE RADIAL
    TYPE (REGIONAL),    DIMENSION(:),   ALLOCATABLE    :: reg
    TYPE(BOUNDARY_CONDITIONS)                          :: bc
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: q
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: wt
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: q_fac
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: inv_q_fac
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: inv_sqrt_q_fac
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: p
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: dp
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: ddp
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: bf
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: dbf
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: ddbf
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: normalization
    REAL(idp)                                          :: end
    REAL(idp)                                          :: sumwt
    INTEGER                                            :: n_reg
    REAL(idp)                                          :: start
    REAL(idp)                                          :: increment
    INTEGER                                            :: first
    INTEGER                                            :: last
    INTEGER                                            :: mmax
    INTEGER                                            :: n_pts
    INTEGER                                            :: n_fun
    INTEGER                                            :: n_phy
    CHARACTER(LEN=16)                                  :: type_quadrature
    INTEGER                                            :: fixed
    INTEGER                                            :: n_fixed
    REAL(idp),          DIMENSION(2)                   :: edge
    INTEGER                                            :: fixed_point
    INTEGER                                            :: n_shl
  END TYPE RADIAL
!
  TYPE PHI
    TYPE (REGIONAL),    DIMENSION(:),    ALLOCATABLE   :: reg
    TYPE(BOUNDARY_CONDITIONS)                          :: bc
    LOGICAL                                            :: fourier
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: q
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: q_fac
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: inv_q_fac
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: inv_sqrt_q_fac
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: p
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: dp
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: ddp
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: bf
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: dbf
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: ddbf
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: normalization
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: wt
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: phi
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: cos_phi
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: sin_phi
    REAL(idp)                                          :: start
    REAL(idp)                                          :: end
    CHARACTER(LEN=16)                                  :: type_quadrature
    REAL(idp)                                          :: increment
    REAL(idp),          DIMENSION(2)                   :: edge
    INTEGER                                            :: fixed_point
    INTEGER                                            :: n_phy
    INTEGER                                            :: n_pts
    INTEGER                                            :: n_reg
    INTEGER                                            :: n_fun
  END TYPE PHI
!
  TYPE ANGULAR                   
    TYPE(LEBEDEV)                                      :: leb_ang
    TYPE(THETA)                                        :: theta_ang
    TYPE(PHI)                                          :: phi_ang
    INTEGER                                            :: n_reg
    INTEGER                                            :: n_phy
    INTEGER                                            :: n_pts
    REAL(idp)                                          :: start
    REAL(idp)                                          :: end
    REAL(idp)                                          :: increment
  END TYPE ANGULAR
!
  TYPE SHELLS
    TYPE(RADIAL)                                       :: rad
    TYPE(ANGULAR)                                      :: ang
    TYPE(LEBEDEV)                                      :: leb_ang
    TYPE(THETA)                                        :: theta_ang
    TYPE(PHI)                                          :: phi_ang
    INTEGER                                            :: n_reg    
    INTEGER                                            :: num
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: q
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: q_fac
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: inv_q_fac
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: inv_sqrt_q_fac
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: p
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: dp
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: ddp
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: bf
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: dbf
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: ddbf
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: normalization
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: wt
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: sin_thet
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: xyz_grid
    REAL(idp),          DIMENSION(:,:), ALLOCATABLE    :: sph_grid
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: weight
    REAL(idp),          DIMENSION(:),   ALLOCATABLE    :: integration_weight
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: yukawa
    REAL(idp),          DIMENSION(2)                   :: edge
    REAL(idp)                                          :: start
    REAL(idp)                                          :: end
    CHARACTER(LEN=16)                                  :: type_quadrature
    REAL(idp)                                          :: increment
    INTEGER                                            :: n_pts
    INTEGER                                            :: n_fixed
    INTEGER                                            :: first
    INTEGER                                            :: last
    INTEGER                                            :: n_fun
    INTEGER                                            :: fixed_point
    INTEGER                                            :: n_2d
    INTEGER                                            :: n_3d
    TYPE (even_matrices)                               :: even
    TYPE (odd_matrices)                                :: odd
  END TYPE SHELLS
!
  TYPE ld06
  END TYPE ld06

  TYPE ld14
  END TYPE ld14

  TYPE ld26
  END TYPE ld26

  TYPE ld38
  END TYPE ld38

  TYPE ld50
  END TYPE ld50

  TYPE ld74
  END TYPE ld74

  TYPE ld86
  END TYPE ld86

  TYPE ld110
  END TYPE ld110

  TYPE ld146
  END TYPE ld146

  TYPE ld170
  END TYPE ld170

  TYPE ld194
  END TYPE ld194

  TYPE ld230
  END TYPE ld230

  TYPE ld266
  END TYPE ld266

  TYPE ld302
  END TYPE ld302

  TYPE ld350
  END TYPE ld350

  TYPE ld434
  END TYPE ld434

  TYPE ld590
  END TYPE ld590

  TYPE ld770
  END TYPE ld770

  TYPE ld974
  END TYPE ld974

  TYPE ld1202
  END TYPE ld1202

  TYPE ld1454
  END TYPE ld1454

  TYPE ld1730
  END TYPE ld1730

  TYPE ld2030
  END TYPE ld2030

  TYPE ld2354
  END TYPE ld2354

  TYPE ld2702
  END TYPE ld2702

  TYPE ld3074
  END TYPE ld3074

  TYPE ld3470
  END TYPE ld3470

  TYPE ld3890
  END TYPE ld3890

  TYPE ld4334
  END TYPE ld4334

  TYPE ld4802
  END TYPE ld4802

  TYPE ld5294
  END TYPE ld5294

  TYPE ld5810
  END TYPE ld5810

  TYPE gl
  END TYPE gl
!
  TYPE RULE     
    INTEGER                                            :: maxang
    Type(gl)                                           :: glrule
    Type(ld06)                                         :: rule_06
    Type(ld14)                                         :: rule_14
    Type(ld26)                                         :: rule_26
    Type(ld38)                                         :: rule_38
    Type(ld50)                                         :: rule_50
    Type(ld74)                                         :: rule_74
    Type(ld86)                                         :: rule_86
    Type(ld110)                                        :: rule_110
    Type(ld146)                                        :: rule_146
    Type(ld170)                                        :: rule_170
    Type(ld194)                                        :: rule_194
    Type(ld230)                                        :: rule_230
    Type(ld266)                                        :: rule_266
    Type(ld302)                                        :: rule_302
    Type(ld350)                                        :: rule_350
    Type(ld434)                                        :: rule_434
    Type(ld590)                                        :: rule_590
    Type(ld770)                                        :: rule_770
    Type(ld974)                                        :: rule_974
    Type(ld1202)                                       :: rule_1202
    Type(ld1454)                                       :: rule_1454
    Type(ld1730)                                       :: rule_1730
    Type(ld2030)                                       :: rule_2030
    Type(ld2354)                                       :: rule_2354
    Type(ld2702)                                       :: rule_2702
    Type(ld3074)                                       :: rule_3074
    Type(ld3470)                                       :: rule_3470
    Type(ld3890)                                       :: rule_3890
    Type(ld4334)                                       :: rule_4334
    Type(ld4802)                                       :: rule_4802
    Type(ld5294)                                       :: rule_5294
    Type(ld5810)                                       :: rule_5810
  END TYPE RULE
!
  TYPE(RULE)                                           :: leb_rule
!
  TYPE(RULE)                                           :: gauss_rule
!
  TYPE WORKING_ARRAYS
  TYPE (FINAL_MATRICES), DIMENSION(:), ALLOCATABLE     :: mat
  REAL(idp),        DIMENSION(:),   ALLOCATABLE        :: eigenvalues
  REAL(idp),        DIMENSION(:),   ALLOCATABLE        :: work
  REAL(idp),        DIMENSION(:),   ALLOCATABLE        :: lower_mat
  REAL(idp),        DIMENSION(:,:), ALLOCATABLE        :: eigenvectors
  REAL(idp),        DIMENSION(2)                       :: edge
  INTEGER                                              :: n_pts
  END TYPE WORKING_ARRAYS
!
  TYPE ATOMS
    TYPE(SHELLS),   DIMENSION(:), ALLOCATABLE          :: shl
    TYPE(REGIONAL), DIMENSION(:), ALLOCATABLE          :: reg
    TYPE(RADIAL)                                       :: rad
    TYPE(ANGULAR)                                      :: ang
    TYPE(LEBEDEV)                                      :: leb_ang
    TYPE(THETA)                                        :: theta_ang
    TYPE(PHI)                                          :: phi_ang
    INTEGER                                            :: n_coord
    INTEGER                                            :: n_shl
    INTEGER                                            :: n_pts
    INTEGER                                            :: shl_no
    INTEGER                                            :: ngrid
    INTEGER                                            :: nwts
    INTEGER                                            :: maxang
    REAL(idp)                                          :: yukawa_integral
    REAL(idp)                                          :: vol
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: a
    REAL(idp)                                          :: eta
    REAL(idp),          DIMENSION(3)                   :: cen
    REAL(idp)                                          :: znuc
    REAL(idp)                                          :: R_max
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: points
    REAL(idp),          DIMENSION(:),    ALLOCATABLE   :: weights
    CHARACTER(LEN=8)                                   :: coord_label
    REAL(isp),          DIMENSION(:),    ALLOCATABLE   :: local_time
    REAL(isp),          DIMENSION(:),    ALLOCATABLE   :: total_time
  END TYPE ATOMS
!
  TYPE(ATOMS),          DIMENSION(:),    ALLOCATABLE   :: atom
  TYPE(ATOMS),          DIMENSION(:),    ALLOCATABLE   :: center
!
  TYPE(WORKING_ARRAYS), DIMENSION(:),    ALLOCATABLE   :: wa
!
END MODULE Grid_Defined_Types




