! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Modification of Recursion Coefficients for Lobatto Quadrature}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck modab.f
  SUBROUTINE modab(a,b,kpts,endpts,n)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: a, b
  INTEGER                                :: kpts
  REAL*8, DIMENSION(2)                   :: endpts
  REAL*8                                 :: gbslve, gam, t1
  IF (kpts == 0) THEN
      RETURN
  ELSE IF (kpts == 1) THEN
!         only a(n) must be changed
      a(n) =gbslve(endpts(1), n, a, b)*b(n-1)**2 + endpts(1)
      RETURN
  ELSE IF(kpts == 2) THEN
!         a(n) and b(n-1) must be recomputed
      gam =gbslve(endpts(1), n, a, b)
      t1 = ((endpts(1) - endpts(2))/(gbslve(endpts(2), n, a, b) - gam))
      b(n-1) =  SQRT(t1)
      a(n) = endpts(1) + gam*t1
      RETURN
  END IF
END SUBROUTINE modab
