!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE r_matrix_global
!deck r_matrix_global.f90
! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{R_MATRIX_GLOBAL: MODULE for R-matrix Code}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!**begin prologue     r_matrix_global
!**date written       030526   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            data entry for dvr r-matrix program.
!**references
!**routines called
!**end prologue       r_matrix_global
  USE io
  USE grid_global
  IMPLICIT NONE
  CHARACTER (LEN=80)                     :: eunit
  LOGICAL, DIMENSION(6)                  :: log_main
  CHARACTER (LEN=80), DIMENSION(6)       :: pr_main
  CHARACTER (LEN=80)                     :: title
  INTEGER                                :: nen, ltop=100           
  REAL*8, DIMENSION(:),   ALLOCATABLE    :: energy
  REAL*8                                 :: rmat, k, rho
  REAL*8, DIMENSION(:),   ALLOCATABLE    :: jbes, djbes,  &
                                            ddjbes, ybes, &
                                            dybes, ddybes
  REAL*8, DIMENSION(:), ALLOCATABLE      :: kmat, phase

  DATA pr_main / 'eigenvalues','eigenvectors',  &
                 'potential','surface-amplitudes', &
                 'r-matrix','none' /
END MODULE r_matrix_global




