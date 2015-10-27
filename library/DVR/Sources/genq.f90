! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Points, Weights and Traansformation Matrix for
!        Generalized Gauss Quadratures}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck @(#)genq.f 1.1 9/9/91
  SUBROUTINE genq( a, b, w, z, fix, edge, muzero, n)
  USE input_output
  IMPLICIT NONE
!        a, b     on input the recursion coefficients; will be destroyed
!     output parameters (both arrays of length n)
!        a        the desired points
!        w        the desired weights
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: a, b, w
  REAL*8, DIMENSION(n,n)                 :: z
  LOGICAL, DIMENSION(2)                  :: fix
  REAL*8, DIMENSION(2)                   :: edge
  REAL*8                                 :: muzero, sumwt
  INTEGER                                :: i, ierr
!     the method used is a ql-type method with origin shifting
  z(:,1)=b
  b(2:n)=z(1:n-1,1)
  w=0.d0
  z=0.d0
  DO  i=1,n
      z(i,i)=1.d0
  END DO
  CALL imtql2 (n, n, a, b, z, ierr)
  w=muzero*z(1,:)*z(1,:)
  sumwt=sum(w,1)
  IF(fix(1)) THEN
     a(1)=edge(1)
  END IF
  IF(fix(2)) THEN
     a(n)=edge(2)
  END IF
  WRITE(iout,1) sumwt
1    FORMAT(/,1X,'sum of the weights = ',e15.8)
END SUBROUTINE genq
