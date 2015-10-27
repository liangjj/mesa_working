! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Final Regional DVR Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck ke_region.f
!***begin prologue     ke_region
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            regional matrix elements of kinetic energy
!                      and potential matrix elements.
!***
!***references
!***routines called
!***end prologue       ke_region
!
  SUBROUTINE ke_region(kmat,nmat,idim)
  USE dvr_global
  USE dvrprop_global
  IMPLICIT NONE
  INTEGER                                :: nmat, idim
  REAL*8, DIMENSION(nmat,nmat)           :: kmat
  INTEGER                                :: start, reg
!
!     calculate the needed matrix elements
!
  start=1
  DO  reg=1,nreg
      CALL fil_reg(kmat(start,start),mat_reg(reg,idim)%ke_mat, &
                   npt(reg),nphy,reg)
      start = start + npt(reg) - 1
  END DO
END SUBROUTINE ke_region



