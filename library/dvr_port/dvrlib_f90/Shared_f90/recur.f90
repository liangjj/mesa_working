! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Recursion coefficients for classical orthogonal Functions}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck recur.f
  SUBROUTINE recur( ckind, n, alpha, beta, kpts, endpts, a, b, muzero )
  USE inout
!           this set of routines computes a and b recursion coefficients
!           for known gauss quadratures

!       input parameters

!        ckind =  legendre: w(x) = 1 on (-1, 1)
!        ckind =  chebyshev-1: chebyshev quadrature of the first kind
!                 w(x) = 1/dsqrt(1 - x*x) on (-1, +1)
!        ckind =  chebyshev-2 quadrature of the second kind
!                 w(x) = dsqrt(1 - x*x) on (-1, 1)
!        ckind =  hermite: w(x) = exp(-x*x) on  (-infinity, +infinity)
!        ckind =  jacobi:  w(x) = (1-x)**alpha * (1+x)**beta
!                 on (-1, 1), alpha, beta .gt. -1.
!                   note= kind=2 and 3 are a special case of this.
!        ckind =  laguerre: w(x) = exp(-x)*x**alpha on
!                 (0, +infinity), alpha .gt. -1

!        n        the number of points used for the quadrature rule
!        alpha    real parameter used only for gauss-jacobi and gauss-
!                 laguerre quadrature (otherwise use 0.).
!        beta     real parameter used only for gauss-jacobi quadrature--
!                 (otherwise use 0.).
!        kpts     (integer) normally 0, unless the left or right end-
!                 point (or both) of the interval is required to be a
!                 node (this is called gauss-radau or gauss-lobatto
!                 quadrature).  then kpts is the number of fixed
!                 endpoints (1 or 2).
!        endpts   real array of length 2.  contains the values of
!                 any fixed endpoints, if kpts = 1 or 2.
!        a        real array of a recursion coefficients
!        b        real array of b recursion coefficients

!     subroutines required

!        gbslve and class are needed

!     ..................................................................
  CHARACTER (LEN=*)          :: ckind
  INTEGER                    :: n, kpts
  REAL*8                     :: alpha, beta
  REAL*8, DIMENSION(2)       :: endpts
  REAL*8, DIMENSION(n)       :: a, b
  REAL*8                     :: muzero
  REAL*8                     :: gbslve, gam, t1
  CHARACTER (LEN=80)         :: title
  a=0.d0
  b=0.d0
  IF(ckind == 'legendre'.or.ckind == 'one') THEN
     kind=1
  ELSE IF(ckind == 'chebyshev-1') THEN
     kind=2
  ELSE IF(ckind == 'chebyshev-2') THEN
     kind=3
  ELSE IF(ckind == 'hermite') THEN
     kind=4
  ELSE IF(ckind == 'jacobi') THEN
     kind=5
  ELSE IF(ckind == 'laguerre') THEN
     kind=6
  ELSE
     write(output,1)
  END IF
!     This routine returns $\alpha$ and $\beta$ for the polynomials.
  CALL class (kind, n, alpha, beta, b, a, muzero)
!     the matrix of coefficients is assumed to be symmetric.
!     make appropriate changes in the lower right 2 by 2 if lobatto rules
!     used.
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
!     note that the indices of the elements of b run from 1 to n-1
!     and thus the value of b(n) is arbitrary.
1 FORMAT('error in quadrature type')
END SUBROUTINE recur


