! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Cylindrical Bloch Operator Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck blcyl.f
!***begin prologue     blcyl
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           bloch operator matrix elements
!***author             schneider, barry (nsf)
!***source
!***purpose            generate bloch operator matrix elements
!***                   in cylindrical coordinates.
!***
!***description
!***references
!***routines called
!***end prologue       blcyl
  SUBROUTINE blcyl(blmat,p,dp,xl,el,xr,er,n,parity)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n,n)                 :: blmat, p, dp
  REAL*8                                 :: xl, el, xr, er
  CHARACTER (LEN=*)                      :: parity
  REAL*8                                 :: two=2.d0
  IF(parity == 'none') THEN
     IF(xl == el) THEN
        blmat(1,:) = blmat(1,:) + xl * p(1,1) *dp(1,:)
     END IF
     IF(xr == er) THEN
        blmat(n,:) = blmat(n,:) - xr * p(n,n) * dp(n,:)
     END IF
  ELSE IF(parity == 'even'.OR.parity == 'odd') THEN
     IF(xl == el) THEN
        blmat(1,:) = blmat(1,:) + two * xl * xl * p(1,1) * dp(1,:)
     END IF
     IF(xr == er) THEN
        blmat(n,:) = blmat(n,:) - two * xr * xr * p(n,n) * dp(n,:)
     END IF
  ELSE
     CALL lnkerr('error in argument type')
  END IF
END SUBROUTINE blcyl
