c \documentclass{article}
c \usepackage{graphicx}
c \setkeys{Gin}{width=\linewidth}
c \title{Recursion coefficients for classical orthogonal Functions}
c \author{Barry I. Schneider}
c \date{}
c \def \<{\langle}
c \def \>{\rangle}
c \begin{document}
c \maketitle    

*deck recur.f
      subroutine recur( ckind, n, alpha, beta, kpts, endpts, 
     1                  a, b, muzero )
c
c           this set of routines computes a and b recursion coefficients
c           for known gauss quadratures
c
c       input parameters
c
c        ckind =  legendre: w(x) = 1 on (-1, 1)
c        ckind =  chebyshev-1: chebyshev quadrature of the first kind
c                 w(x) = 1/dsqrt(1 - x*x) on (-1, +1)
c        ckind =  chebyshev-2 quadrature of the second kind
c                 w(x) = dsqrt(1 - x*x) on (-1, 1)
c        ckind =  hermite: w(x) = exp(-x*x) on  (-infinity, +infinity)
c        ckind =  jacobi:  w(x) = (1-x)**alpha * (1+x)**beta
c                 on (-1, 1), alpha, beta .gt. -1.
c                   note= kind=2 and 3 are a special case of this.
c        ckind =  laguerre: w(x) = exp(-x)*x**alpha on 
c                 (0, +infinity), alpha .gt. -1
c
c        n        the number of points used for the quadrature rule
c        alpha    real parameter used only for gauss-jacobi and gauss-
c                 laguerre quadrature (otherwise use 0.).
c        beta     real parameter used only for gauss-jacobi quadrature--
c                 (otherwise use 0.).
c        kpts     (integer) normally 0, unless the left or right end-
c                 point (or both) of the interval is required to be a
c                 node (this is called gauss-radau or gauss-lobatto
c                 quadrature).  then kpts is the number of fixed
c                 endpoints (1 or 2).
c        endpts   real array of length 2.  contains the values of
c                 any fixed endpoints, if kpts = 1 or 2.
c        a        real array of a recursion coefficients
c        b        real array of b recursion coefficients
c
c     subroutines required
c
c        gbslve and class are needed
c
c     ..................................................................
c
      implicit integer(a-z)
      real*8 alpha, beta, a, b, endpts, muzero
      real*8 gbslve, gam, t1
      character *(*) ckind
      character*80 title
      dimension  a(n), b(n), endpts(2)
      common/io/inp,iout
      call rzero(a,n)
      call rzero(b,n)
      if(ckind.eq.'legendre') then
          kind=1
      elseif(ckind.eq.'chebyshev-1') then
          kind=2
      elseif(ckind.eq.'chebyshev-2') then
          kind=3
      elseif(ckind.eq.'hermite') then
          kind=4
      elseif(ckind.eq.'jacobi') then
          kind=5
      elseif(ckind.eq.'laguerre') then
          kind=6
      else
          call lnkerr('error in quadrature type')
      endif

c     This routine returns $\alpha$ and $\beta$ for the polynomials.

      call class (kind, n, alpha, beta, b, a, muzero)

c     the matrix of coefficients is assumed to be symmetric.
c     make appropriate changes in the lower right 2 by 2 if lobatto rules
c     used.

      if (kpts.eq.0) then
          return
      elseif (kpts.eq.1) then

c         only a(n) must be changed

          a(n) =gbslve(endpts(1), n, a, b)*b(n-1)**2 + endpts(1)
          return
      elseif(kpts.eq.2) then

c         a(n) and b(n-1) must be recomputed
 
          gam =gbslve(endpts(1), n, a, b)
          t1 = ((endpts(1) - endpts(2))/(gbslve(endpts(2), n, a, b) 
     1                     - gam))
          b(n-1) =  sqrt(t1)
          a(n) = endpts(1) + gam*t1
          return
      endif

c     note that the indices of the elements of b run from 1 to n-1
c     and thus the value of b(n) is arbitrary.

      return
      end


