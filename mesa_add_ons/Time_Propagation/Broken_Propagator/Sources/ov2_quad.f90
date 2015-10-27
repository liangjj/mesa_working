! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Subroutine for Gaussian Wavepacket Normalization}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck ov2_quad
!**begin prologue     ov2_quad
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate normalization integrals for 
!***                  gaussian wavepacket
!**references
!**routines called
!**end prologue       ov2_quad
  SUBROUTINE ov2_quad(ov,pt,wt,x_0,alpha,sigma,n)
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: pt, wt
  REAL*8                                 :: ov, x_0, alpha, sigma
  INTEGER                                :: i 
  ov=0.d0
  DO  i=1,n
      ov = ov + wt(i) * pt(i) * pt(i) * EXP(- alpha * ( pt(i) - x_0 ) *  &
        ( pt(i) - x_0 ) / ( sigma * sigma ) )
  END DO
END SUBROUTINE ov2_quad







