!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE prop_global
!deck prop_global.f
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
!**begin prologue     prop_global
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            global arrays for time propagator.
!**references
!**routines called
!**end prologue       prop_global
  USE grid_global
  USE potential
  IMPLICIT NONE
  CHARACTER (LEN=24)                     :: typ_pak, plot_name
  CHARACTER (LEN=80)                     :: algorithm
  REAL*8                                 :: cnverg, thresh
  REAL*4, DIMENSION(20)                  :: delta, time
  INTEGER, DIMENSION(6)                  :: iplot
  INTEGER                                :: ntreg
  INTEGER                                :: ndiff
  INTEGER                                :: n3d
  INTEGER                                :: nchan
  INTEGER                                :: vec_size
  INTEGER                                :: maxit, maxvec, ntrial
  INTEGER                                :: plot_step
  INTEGER                                :: state
  INTEGER, DIMENSION(3)                  :: pun
  LOGICAL                                :: imtime, exact, restart
  LOGICAL                                :: space
  LOGICAL                                :: plot, proj, prnton
  CHARACTER (LEN=80)                     :: i0stat, title
  REAL*8,   DIMENSION(:),                                              &
            ALLOCATABLE                  :: v_tot
  REAL*8,   DIMENSION(:),                                              &
            ALLOCATABLE                  :: tim_pts, tedge
  COMPLEX*16, DIMENSION(:),                                            &
            ALLOCATABLE                  :: auto_corr
  COMPLEX*16                             :: eye
  REAL*8                                 :: vt, energy
  REAL*8                                 :: t0, t1, deltat
  LOGICAL                                :: nochk
  LOGICAL, DIMENSION(9)                  :: log_prp
  LOGICAL, DIMENSION(9)                  :: log_main
  DATA iplot / 95, 96, 97, 98, 99, 100 /
  DATA eye / (0.d0,1.d0) /
  DATA plot_name / 'wave_function' /
END MODULE prop_global




