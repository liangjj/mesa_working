!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE dvrprop_global_it
!deck dvrprop_global_it.f
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
!**begin prologue     dvrprop_global_it
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            global arrays for FEDVR time propagator.
!**references
!**routines called
!**end prologue       dvrprop_global_it
  USE input_output
  USE prop_global
  USE prop_prnt
  IMPLICIT NONE
  CHARACTER(LEN=8)                               :: key
  TYPE dvr_mat_d
    REAL*8,         DIMENSION(:,:),                                          &
                    POINTER                      :: ke_mat_d, eigvec_mat_d
    REAL*8,         DIMENSION(:),                                            &
                    POINTER                      :: pt_d, eigval_mat_d
    REAL*8,         DIMENSION(:,:,:),                                        &
                    POINTER                      :: exp_t_mat
  END TYPE dvr_mat_d
  TYPE (dvr_mat_d), DIMENSION(:,:),                                          &
                    ALLOCATABLE                  :: mat_reg_d
  REAL*8,           DIMENSION(:),                                            &
                    ALLOCATABLE                  :: psi_1d, v_scr_1d 
  REAL*8,           DIMENSION(:,:),                                          &
                    ALLOCATABLE                  :: psi_2d, v_scr_2d 
  REAL*8,           DIMENSION(:,:,:)  ,                                      &
                    ALLOCATABLE                  :: psi_3d, v_scr_3d
  REAL*8,           DIMENSION(:),                                            &
                    ALLOCATABLE                  :: exp_diag
  REAL*8,           DIMENSION(:),                                            &
                    ALLOCATABLE                  :: scr_d
  REAL*8,           DIMENSION(:,:),                                          &
                    ALLOCATABLE                  :: exp_tmp_d
  REAL*8                                         :: t_init
  TYPE c_temp
    REAL*8,         DIMENSION(:),                                            &
                    POINTER                      :: c, phi
    INTEGER,        DIMENSION(:),                                            &
                    POINTER                      :: l
  END TYPE c_temp
  TYPE(c_temp),     DIMENSION(:),                                            &
                    ALLOCATABLE                  :: c_loc
  TYPE z_temp
    REAL*8,         DIMENSION(:),                                            &
                    POINTER                      :: z
  END TYPE z_temp
  TYPE(z_temp),     DIMENSION(:),                                            &
                    ALLOCATABLE                  :: zloc
  INTEGER                                        :: maxmem_reg, prop_order
  INTEGER                                        :: n_prop, prop_point
  INTEGER                                        :: maxdim  
  REAL*8                                         :: p_fac, p_loc, tau_loc
  INTEGER                                        :: nvec  
  CHARACTER (LEN=16)                             :: e_method
  CHARACTER (LEN=24)                             :: diag_mod
  REAL*8,          DIMENSION(:),                                             &
                   ALLOCATABLE                   :: f_1, fac, v_1
  REAL*8,          DIMENSION(:,:),                                           &
                   ALLOCATABLE                   :: f_2
  REAL*8,          DIMENSION(:,:,:),                                         &
                   ALLOCATABLE                   :: f_3
  REAL*8                                         :: con
END MODULE dvrprop_global_it




