!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE dvrprop_global_rt
!deck dvrprop_global_rt.f
! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Prop_GLOBAL: MODULE for Arnoldi Propagation}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!**begin prologue     dvrprop_global_rt
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            global arrays for FEDVR time propagator.
!**references
!**routines called
!**end prologue       dvrprop_global_rt
  USE input_output
  USE prop_global
  USE prop_prnt
  IMPLICIT NONE
  CHARACTER(LEN=8)                                :: key
  CHARACTER (LEN=24)                              :: diag_mod
  LOGICAL                                         :: absorb
  INTEGER, DIMENSION(4)                           :: n_reg_real, n_reg_absorb
  INTEGER                                         :: starting_reg, ending_reg, n_reg
  TYPE dvr_mat_d
    REAL*8,          DIMENSION(:,:),                                         &
                     POINTER                      :: ke_mat_d, eigvec_mat_d
    REAL*8,          DIMENSION(:),                                           &
                     POINTER                      :: pt_d, eigval_mat_d
    REAL*8,          DIMENSION(:,:,:),                                       &
                     POINTER                      :: cosine_t_mat, sine_t_mat
    COMPLEX*16,      DIMENSION(:,:),                                         &
                     POINTER                      :: exp_t_mat
  END TYPE dvr_mat_d
  TYPE dvr_mat_z
    COMPLEX*16,      DIMENSION(:,:),                                         &
                     POINTER                      :: ke_mat_z, eigvec_mat_z_l
    COMPLEX*16,      DIMENSION(:,:),                                         &
                     POINTER                      :: eigvec_mat_z_r
    COMPLEX*16,      DIMENSION(:),                                           &
                     POINTER                      :: eigval_mat_z
    COMPLEX*16,      DIMENSION(:),                                           &
                     POINTER                      :: v_add_z
  END TYPE dvr_mat_z
  TYPE (dvr_mat_d),  DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: mat_reg_d
  TYPE (dvr_mat_z),  DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: mat_reg_z
  COMPLEX*16,        DIMENSION(:),                                           &
                     ALLOCATABLE                  :: psi_z, v_scr_z 
  REAL*8,            DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: psi_1d, v_scr_1d 
  REAL*8,            DIMENSION(:,:,:),                                       &
                     ALLOCATABLE                  :: psi_2d, v_scr_2d 
  REAL*8,            DIMENSION(:,:,:,:),                                     &
                     ALLOCATABLE                  :: psi_3d, v_scr_3d 
  REAL*8,            DIMENSION(:),                                           &
                     ALLOCATABLE                  :: sin_diag, cos_diag, exp_diag
  COMPLEX*16,        DIMENSION(:),                                           &
                     ALLOCATABLE                  :: exp_diag_z
  REAL*8,            DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: si_d, ci_d
  COMPLEX*16,        DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: si_z, ci_z
  REAL*8,            DIMENSION(:),                                           &
                     ALLOCATABLE                  :: scr_d
  REAL*8,            DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: exp_tmp_d
  COMPLEX*16,        DIMENSION(:),                                           &
                     ALLOCATABLE                  :: scr_z
  REAL*8                                          :: t_init
  TYPE c_temp
    REAL*8,          DIMENSION(:),                                           &
                     POINTER                      :: c, phi
    INTEGER,         DIMENSION(:),                                           &
                     POINTER                      :: l
  END TYPE c_temp
  TYPE(c_temp),      DIMENSION(:),                                           &
                     ALLOCATABLE                  :: c_loc
  TYPE z_temp
    COMPLEX*16,      DIMENSION(:),                                           &
                     POINTER                      :: z
  END TYPE z_temp
  TYPE(z_temp),      DIMENSION(:),                                           &
                     ALLOCATABLE                  :: zloc
  INTEGER                                         :: maxmem_reg, prop_order
  INTEGER                                         :: n_prop, prop_point
  INTEGER                                         :: maxdim  
  INTEGER                                         :: nr_absorb  
  REAL*8                                          :: p_fac, p_loc, tau_loc
  INTEGER                                         :: nvec  
  REAL*8,            DIMENSION(:),                                           &
                     ALLOCATABLE                  :: f_1, fac, v_1
  REAL*8,            DIMENSION(:,:),                                         &
                     ALLOCATABLE                  :: f_2
  REAL*8,            DIMENSION(:,:,:),                                       &
                     ALLOCATABLE                  :: f_3
  REAL*8                                          :: con
END MODULE dvrprop_global_rt




