! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Cartesian Bloch Operator Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck blxyz.f
!***begin prologue     blxyz
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           bloch matrix elements
!***author             schneider, barry (nsf)
!***source
!***purpose            generate bloch matrix elements
!***                   in cartesian coordinates.
!***
!***description
!***references

!***routines called
!***end prologue       blxyz
  SUBROUTINE blxyz(blmat,p,dp,xl,el,xr,er,n)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n    
  REAL*8, DIMENSION(n,n)                 :: blmat, p, dp
  REAL*8                                 :: xl, el, xr, er
  IF(xl == el) THEN
     blmat(1,:) = blmat(1,:) + p(1,1)*dp(1,:)
  END IF
  IF(xr == er) THEN
     blmat(n,:) = blmat(n,:) - p(n,n)*dp(n,:)
  END IF
END SUBROUTINE blxyz
