! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Cartesian Kinetic Energy Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck kinxyz.f
!***begin prologue     kinxyz
!***date written       022202   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           kinetic energy
!***author             schneider, barry (nsf)
!***source
!***purpose            generate kinetic energy matrix elements
!***                   in cartesian coordinates.
!***
!***description
!***references
!***routines called
!***end prologue       kinxyz
  SUBROUTINE kinxyz(kmat,p,ddp,wt,n)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n    
  REAL*8, DIMENSION(n,n)                 :: kmat, p, ddp
  REAL*8, DIMENSION(n)                   :: wt
  INTEGER                                :: i, j   
  DO  i=1,n
      kmat(i,:) = kmat(i,:) + p(i,i)*ddp(i,:)*wt(i)
  END DO
END SUBROUTINE kinxyz
