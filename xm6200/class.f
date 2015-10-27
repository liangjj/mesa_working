*deck %W%  %G%
      subroutine class(kind, n, alpha, beta, b, a, muzero)
c
c           this procedure supplies the coefficients a(j), b(j) of the
c        recurrence relation
c
c             b p (x) = (x - a ) p   (x) - b   p   (x)
c              j j            j   j-1       j-1 j-2
c
c        for the various classical (normalized) orthogonal polynomials,
c        and the zero-th moment
c
c             muzero = integral w(x) dx
c
c        of the given polynomial   weight function w(x).  since the
c        polynomials are orthonormalized, the tridiagonal matrix is
c        guaranteed to be symmetric.
c
c           the input parameter alpha is used only for laguerre and
c        jacobi polynomials, and the parameter beta is used only for
c        jacobi polynomials.  the laguerre and jacobi polynomials
c        require the gamma function.
c
c     ..................................................................
c
      implicit real*8 (a-h,o-z)
      dimension  a(n),b(n)
      common/io/inp,iout
      real*8  muzero
      data pi / 3.141592653589793d0  /
c
      nm1 = n - 1
      go to (10, 20, 30, 40, 50, 60), kind
c
c              kind = 1=  legendre polynomials p(x)
c              on (-1, +1), w(x) = 1.
c
   10 muzero = 2.0d0
      do 11 i = 1, nm1
         a(i) = 0.0d0
         abi = i
   11    b(i) = abi/ sqrt(4.d0*abi*abi - 1.0d0  )
      a(n) = 0.0d0
      return
c
c              kind = 2=  chebyshev polynomials of the first kind t(x)
c              on (-1, +1), w(x) = 1 / sqrt(1 - x*x)
c
   20 muzero = pi
      do 21 i = 1, nm1
         a(i) = 0.0d0
   21    b(i) = 0.5d0
      b(1) =  sqrt(0.5d0  )
      a(n) = 0.0d0
      return
c
c              kind = 3=  chebyshev polynomials of the second kind u(x)
c              on (-1, +1), w(x) = sqrt(1 - x*x)
c
   30 muzero = pi/2.0d0
      do 31 i = 1, nm1
         a(i) = 0.0d0
   31    b(i) = 0.5d0
      a(n) = 0.0d0
      return
c
c              kind = 4=  hermite polynomials h(x) on (-infinity,
c              +infinity), w(x) = exp(-x**2)
c
   40 muzero =  sqrt(pi)
      do 41 i = 1, nm1
         a(i) = 0.0d0
   41    b(i) =  sqrt(i/2.0d0  )
      a(n) = 0.0d0
      return
c
c              kind = 5=  jacobi polynomials p(alpha, beta)(x) on
c              (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta, alpha and
c              beta greater than -1
c
   50 ab = alpha + beta
      abi = 2.0d0   + ab
      muzero = 2.0d0   ** (ab + 1.0d0  ) * gamfun(alpha + 1.0d0  ) * gam
     vfun(
     x beta + 1.0d0  ) / gamfun(abi)
      a(1) = (beta - alpha)/abi
      b(1) =  sqrt(4.0d0  *(1.0d0  + alpha)*(1.0d0   + beta)/((abi + 1.
     v0d0  )*
     1  abi*abi))
      a2b2 = beta*beta - alpha*alpha
      do 51 i = 2, nm1
         abi = 2.0d0  *i + ab
         a(i) = a2b2/((abi - 2.0d0  )*abi)
   51    b(i) =  sqrt (4.0d0  *i*(i + alpha)*(i + beta)*(i + ab)/
     1   ((abi*abi - 1)*abi*abi))
      abi = 2.0d0  *n + ab
      a(n) = a2b2/((abi - 2.0d0  )*abi)
      return
c
c              kind = 6=  laguerre polynomials l(alpha)(x) on
c              (0, +infinity), w(x) = exp(-x) * x**alpha, alpha greater
c              than -1.
c
   60 muzero = gamfun(alpha + 1.0d0  )
      do 61 i = 1, nm1
         a(i) = 2.0d0  *i - 1.0d0   + alpha
   61    b(i) =  sqrt(i*(i + alpha))
      a(n) = 2.0d0  *n - 1 + alpha
      return
      end
