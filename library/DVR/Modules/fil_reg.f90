! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Fill Regional Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck fil_reg.f
!***begin prologue     fil_reg
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           coordinate eigenfunctions
!***author             schneider, barry (nsf)
!***source
!***purpose            Fill regional matrices from global ones.
!***
!***description
!***references
!***routines called
!***end prologue       fil_reg
!\begin{eqnarray}
!\end{eqnarray}
  SUBROUTINE fil_reg(kmat,tr,nr,n,reg)
  USE dvr_global,   ONLY  : iout, prn
  IMPLICIT NONE
  INTEGER                                :: nr, n, reg
  REAL*8                                 :: kmat(n,*) 
  REAL*8, DIMENSION(nr,nr)               :: tr
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  INTEGER                                :: i
!
  tr = kmat(1:nr,1:nr)
  DO i=1,nr
   tr(i,i) = 0.d0
  END DO
  IF(prn(12)) then
     title='Kinetic Energy Matrix Minus Diagonals for '// 'region = '//itoc(reg)
     CALL prntrm(title,tr,nr,nr,nr,nr,iout)
  END IF
END SUBROUTINE fil_reg
