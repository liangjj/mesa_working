!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE dvrprop_global
!deck dvrprop_global.f
!**begin prologue     dvrprop_global
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            global arrays for FEDVR time propagator.
!**references
!**routines called
!**end prologue       dvrprop_global
  USE input_output
  USE prop_global
  USE prop_prnt
  IMPLICIT NONE
  CHARACTER(LEN=8)                                :: key
  CHARACTER (LEN=24)                              :: diag_mod
  LOGICAL                                         :: absorb, no_cc, no_pot
  INTEGER, DIMENSION(4)                           :: n_reg_real, n_reg_absorb
  INTEGER                                         :: starting_reg, ending_reg, n_reg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE dvr_mat
    REAL*8,          DIMENSION(:,:),                                         &
                     POINTER                      :: ke_mat_d
    REAL*8,          DIMENSION(:),                                           &
                     POINTER                      :: pt_d
    REAL*8,          DIMENSION(:,:),                                         &
                     POINTER                      :: eigvec_mat_d
    REAL*8,          DIMENSION(:),                                           &
                     POINTER                      :: eigval_mat_d
    REAL*8,          DIMENSION(:,:,:),                                       &
                     POINTER                      :: exp_t_mat_d
    COMPLEX*16,      DIMENSION(:,:,:),                                       &
                     POINTER                      :: exp_t_mat_z
    REAL*8,          DIMENSION(:,:,:),                                       &
                     POINTER                      :: cos_t_mat
    REAL*8,          DIMENSION(:,:,:),                                       &
                     POINTER                      :: sin_t_mat
    COMPLEX*16,      DIMENSION(:,:),                                         &
                     POINTER                      :: ke_mat_z
    COMPLEX*16,      DIMENSION(:),                                           &
                     POINTER                      :: eigval_mat_z
    COMPLEX*16,      DIMENSION(:,:),                                         &
                     POINTER                      :: eigvec_mat_z_l
    COMPLEX*16,      DIMENSION(:,:),                                         &
                     POINTER                      :: eigvec_mat_z_r
    COMPLEX*16,      DIMENSION(:),                                           &
                     POINTER                      :: v_add_z
  END TYPE dvr_mat

  TYPE (dvr_mat),    DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: mat_reg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL*8,            DIMENSION(:),                                           &
                     ALLOCATABLE                  :: psi_d, v_scr_d 
  REAL*8,            DIMENSION(:),                                           &
                     ALLOCATABLE                  :: psi_1d_d, v_scr_1d_d 
  REAL*8,            DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: psi_2d_d, v_scr_2d_d 
  REAL*8,            DIMENSION(:,:,:),                                       &
                     ALLOCATABLE                  :: psi_3d_d, v_scr_3d_d 
  COMPLEX*16,        DIMENSION(:),                                           &
                     ALLOCATABLE                  :: psi_z, v_scr_z 
  REAL*8,            DIMENSION(:),                                           &
                     ALLOCATABLE                  :: psi_r
  REAL*8,            DIMENSION(:),                                           &
                     ALLOCATABLE                  :: psi_i
  COMPLEX*16,        DIMENSION(:),                                           &
                     ALLOCATABLE                  :: psi_1d_z, v_scr_1d_z 
  COMPLEX*16,        DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: psi_2d_z, v_scr_2d_z 
  COMPLEX*16,        DIMENSION(:,:,:),                                       &
                     ALLOCATABLE                  :: psi_3d_z, v_scr_3d_z 
  REAL*8,            DIMENSION(:),                                           &
                     ALLOCATABLE                  :: exp_diag_d
  COMPLEX*16,        DIMENSION(:),                                           &
                     ALLOCATABLE                  :: exp_diag_z
  REAL*8,            DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: exp_d
  COMPLEX*16,        DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: exp_z
  REAL*8,            DIMENSION(:),                                           &
                     ALLOCATABLE                  :: scr_d
  COMPLEX*16,        DIMENSION(:),                                           &
                     ALLOCATABLE                  :: scr_z
  REAL*8,            DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: exp_tmp_d
  REAL*8,            DIMENSION(:),                                           &
                     ALLOCATABLE                  :: etemp
  INTEGER,           DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: itemp
  REAL*8                                          :: t_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE c_temp
    REAL*8,          DIMENSION(:),                                           &
                     POINTER                      :: c, phi
    INTEGER,         DIMENSION(:),                                           &
                     POINTER                      :: l
  END TYPE c_temp
  TYPE(c_temp),      DIMENSION(:),                                           &
                     ALLOCATABLE                  :: c_loc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE z_temp
    REAL*8,          DIMENSION(:),                                           &
                     POINTER                      :: z_d
    COMPLEX*16,      DIMENSION(:),                                           &
                     POINTER                      :: z_z
  END TYPE z_temp
  TYPE(z_temp),      DIMENSION(:),                                           &
                     ALLOCATABLE                  :: zloc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER                                         :: maxmem_reg, prop_order
  INTEGER                                         :: n_prop, prop_point
  INTEGER                                         :: maxdim  
  INTEGER                                         :: nr_absorb  
  REAL*8                                          :: p_fac, p_loc, tau_loc
  INTEGER                                         :: nvec  
  REAL*8                                          :: con
  CHARACTER (LEN=16)                              :: e_method
  CHARACTER (LEN=16)                              :: type_calculation
  REAL*8                                          :: eig_old, eig_tts
  REAL*8                                          :: eig_tst
  REAL*8                                          :: e_cur
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE dvrprop_global




