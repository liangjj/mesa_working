      program gauss
c sample driver to compute Gauss-Lobatto quadrature points and weights
c n is the order of the quadrature
      implicit real*8 (a-h,o-z)
      real*8 b(20),t(20),w(20),endpts(2)
      kind=1
      z=gamfun(3.d0)
      write(6,*)z
      n=20
      kpts=2
      endpts(1)=-1.
      endpts(2)=+1.
      call gaussq(kind,n,alpha,beta,kpts,endpts,b,t,w)
      write(6,100)(t(i),w(i),i=1,n)
 100  format(2e16.8)
      stop
      end
      subroutine gaussq(kind, n, alpha, beta, kpts, endpts, b, t, w)
c
c           this set of routines computes the nodes x(i) and weights
c        c(i) for gaussian-type quadrature rules with pre-assigned
c        nodes.  these are used when one wishes to approximate
c
c                 integral (from a to b)  f(x) w(x) dx
c
c                              n
c        by                   sum c  f(x )
c                             i=1  i    i
c
c        here w(x) is one of six possible non-negative weight
c        functions (listed below), and f(x) is the
c        function to be integrated.  gaussian quadrature is particularly
c        useful on infinite intervals (with appropriate weight
c        functions), since then other techniques often fail.
c
c           associated with each weight function w(x) is a set of
c        orthogonal polynomials.  the nodes x(i) are just the zeroes
c        of the proper n-th degree polynomial.
c
c     input parameters
c
c        kind     an integer between 0 and 6 giving the type of
c                 quadrature rule
c
c        kind = 0=  simpson's rule w(x) = 1 on (-1, 1) n must be odd.
c        kind = 1=  legendre quadrature, w(x) = 1 on (-1, 1)
c        kind = 2=  chebyshev quadrature of the first kind
c                   w(x) = 1/sqrt(1 - x*x) on (-1, +1)
c        kind = 3=  chebyshev quadrature of the second kind
c                   w(x) = sqrt(1 - x*x) on (-1, 1)
c        kind = 4=  hermite quadrature, w(x) = exp(-x*x) on
c                   (-infinity, +infinity)
c        kind = 5=  jacobi quadrature, w(x) = (1-x)**alpha * (1+x)**
c                   beta on (-1, 1), alpha, beta .gt. -1.
c                   note= kind=2 and 3 are a special case of this.
c        kind = 6=  generalized laguerre quadrature, w(x) = exp(-x)*
c                   x**alpha on (0, +infinity), alpha .gt. -1
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
c        b        real scratch array of length n
c
c     output parameters (both arrays of length n)
c
c        t        will contain the desired nodes x(1),,,x(n)
c        w        will contain the desired weights c(1),,,c(n)
c
c     subroutines required
c
c        gbslve, class, and gbtql2 are provided. underflow may sometimes
c        occur, but it is harmless if the underflow interrupts are
c        turned off as they are on this machine.
c
c     accuracy
c
c        the routine was tested up to n = 512 for legendre quadrature,
c        up to n = 136 for hermite, up to n = 68 for laguerre, and up
c        to n = 10 or 20 in other cases.  in all but two instances,
c        comparison with tables in ref. 3 showed 12 or more significant
c        digits of accuracy.  the two exceptions were the weights for
c        hermite and laguerre quadrature, where underflow caused some
c        very small weights to be set to zero.  this is, of course,
c        completely harmless.
c
c     method
c
c           the coefficients of the three-term recurrence relation
c        for the corresponding set of orthogonal polynomials are
c        used to form a symmetric tridiagonal matrix, whose
c        eigenvalues (determined by the implicit ql-method with
c        shifts) are just the desired nodes.  the first components of
c        the orthonormalized eigenvectors, when properly scaled,
c        yield the weights.  this technique is much faster than using a
c        root-finder to locate the zeroes of the orthogonal polynomial.
c        for further details, see ref. 1.  ref. 2 contains details of
c        gauss-radau and gauss-lobatto quadrature only.
c
c     references
c
c        1.  golub, g. h., and welsch, j. h.,  calculation of gaussian
c            quadrature rules,  mathematics of computation 23 (april,
c            1969), pp. 221-230.
c        2.  golub, g. h.,  some modified matrix eigenvalue problems,
c            siam review 15 (april, 1973), pp. 318-334 (section 7).
c        3.  stroud and secrest, gaussian quadrature formulas, prentice-
c            hall, englewood cliffs, n.j., 1966.
c
c     ..................................................................
c
      implicit real*8 (a-h,o-z)
      real*8  muzero
      dimension  b(n),t(n),w(n),endpts(2)
      if(kind.eq.0) then
       if(2*(n/2).eq.n) then
        write(6,800) n
800     format(" n must be odd for simpson's rule ",i5)
        stop
      endif
        if(n.le.1) then
        t(1) = 0.
        w(1) = 2.0
        return
      endif
       h = 2.0/(n-1)
       t(1) = -1.0
       t(n) = 1.0
       w(1) = h/3.0
       w(n) = h/3.0
       nm1 = n-1
       do 801 i=2,nm1
       t(i) = t(i-1) + h
       w(i) = 4.0 - 2.0*(i-2*(i/2))
801    w(i) = w(i)*h/3.0
       return
      endif
c
      call class (kind, n, alpha, beta, b, t, muzero)
c
c           the matrix of coefficients is assumed to be symmetric.
c           the array t contains the diagonal elements, the array
c           b the off-diagonal elements.
c           make appropriate changes in the lower right 2 by 2
c           submatrix.
c
      if (kpts.eq.0)  go to 100
      if (kpts.eq.2)  go to  50
c
c           if kpts=1, only t(n) must be changed
c
      t(n) =gbslve(endpts(1), n, t, b)*b(n-1)**2 + endpts(1)
      go to 100
c
c           if kpts=2, t(n) and b(n-1) must be recomputed
c
   50 gam =gbslve(endpts(1), n, t, b)
      t1 = ((endpts(1) - endpts(2))/(gbslve(endpts(2), n, t, b) - gam))
      b(n-1) =  sqrt(t1)
      t(n) = endpts(1) + gam*t1
c
c           note that the indices of the elements of b run from 1 to n-1
c           and thus the value of b(n) is arbitrary.
c           now compute the eigenvalues of the symmetric tridiagonal
c           matrix, which has been modified as necessary.
c           the method used is a ql-type method with origin shifting
c
  100 w(1) = 1.0d0
      do 105 i = 2, n
  105    w(i) = 0.0d0
c
      call gbtql2 (n, t, b, w, ierr)
      do 110 i = 1, n
  110    w(i) = muzero * w(i) * w(i)
c
      return
      end
c
c
c
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
   11    b(i) = abi/ sqrt(4*abi*abi - 1.0d0  )
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
c
c
c
      function gbslve(shift, n, a, b)
c
c       this procedure performs elimination to solve for the
c       n-th component of the solution delta to the equation
c
c             (jn - shift*identity) * delta  = en,
c
c       where en is the vector of all zeroes except for 1 in
c       the n-th position.
c
c       the matrix jn is symmetric tridiagonal, with diagonal
c       elements a(i), off-diagonal elements b(i).  this equation
c       must be solved to obtain the appropriate changes in the lower
c       2 by 2 submatrix of coefficients for orthogonal polynomials.
c
c
      implicit real*8 (a-h,o-z)
      dimension  a(n),b(n)
c
      alpha = a(1) - shift
      nm1 = n - 1
      do 10 i = 2, nm1
   10    alpha = a(i) - shift - b(i-1)**2/alpha
      gbslve = 1.0d0  /alpha
      return
      end
c     ------------------------------------------------------------------
c
      subroutine gbtql2(n, d, e, z, ierr)
c
c     this subroutine is a translation of the algol procedure imtql2,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues and first components of the
c     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
c     method, and is adapted from the eispak routine imtql2
c
c     on input=
c
c        n is the order of the matrix;
c
c        d contains the diagonal elements of the input matrix;
c
c        e contains the subdiagonal elements of the input matrix
c          in its first n-1 positions.  e(n) is arbitrary;
c
c        z contains the first row of the identity matrix.
c
c      on output=
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1, 2, ..., ierr-1;
c
c        e has been destroyed;
c
c        z contains the first components of the orthonormal eigenvectors
c          of the symmetric tridiagonal matrix.  if an error exit is
c          made, z contains the eigenvectors associated with the stored
c          eigenvalues;
c
c        ierr is set to
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     ------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      integer i, j, k, l, m, n, ii, mml, ierr
      real*8  machep
      dimension  d(n),e(n),z(n)
c
c     ========== machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c                machep = 16.0d0**(-13) for long form arithmetic
c                on s360 ==========
       machep=1.0e-14
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      e(n) = 0.0d0
      do 240 l = 1, n
         j = 0
c     ========== look for small sub-diagonal element ==========
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            if ( abs(e(m)) .le. machep * ( abs(d(m)) +  abs(d(m+1))))
     x         go to 120
  110    continue
c
  120    p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
c     ========== form shift ==========
         g = (d(l+1) - p) / (2.0d0   * e(l))
         r =  sqrt(g*g+1.0d0  )
         g = d(m) - p + e(l) / (g +  sign(r, g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c     ========== for i=m-1 step -1 until l do -- ==========
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if ( abs(f) .lt.  abs(g)) go to 150
            c = g / f
            r =  sqrt(c*c+1.0d0  )
            e(i+1) = f * r
            s = 1.0d0   / r
            c = c * s
            go to 160
  150       s = f / g
            r =  sqrt(s*s+1.0d0  )
            e(i+1) = g * r
            c = 1.0d0   / r
            s = s * c
  160       g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0   * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
c     ========== form first component of vector ==========
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
            z(i) = c * z(i) - s * f
c
  200    continue
c
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
  240 continue
c     ========== order eigenvalues and eigenvectors ==========
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         p = z(i)
         z(i) = z(k)
         z(k) = p
c
  300 continue
c
      go to 1001
c     ========== set error -- no convergence to an
c                eigenvalue after 30 iterations ==========
 1000 ierr = l
 1001 return
c     ========== last card of gbtql2 ==========
      end
      double precision function  gamfun(z)
c  this is a procedure that evaluates gamma(z) for
c     0 lt z le 3 to 16 significant figures
c    it is based on a chebyshev-type polynomial
c   approximation given in h. werner and r. collinge, math. comput.
c    15 (1961), pp. 195-97.
c   approximations to the gamma function, accurate up to 18 significant
c   digits, may be found in the paper quoted above
c
c
c
      implicit real*8 (a-h,o-z)
      dimension  a(18)
c
       a(1)=1.0d0
       a(2)=.4227843350984678d0
       a(3)=.4118403304263672d0
      a(4)=.0815769192502609d0
      a(5)=.0742490106800904d0
      a(6)=-.0002669810333484d0
      a(7)=.0111540360240344d0
      a(8)=-.0028525821446197d0
      a(9)=.0021036287024598d0
      a(10)=-.0009184843690991d0
      a(11)=.0004874227944768d0
      a(12)=-.0002347204018919d0
      a(13)=.0001115339519666d0
      a(14)=-.0000478747983834d0
      a(15)=.0000175102727179d0
      a(16)=-.0000049203750904d0
      a(17)=.0000009199156407d0
      a(18)=-.0000000839940496d0
c
c
c
      if(z.le.1.0d0  ) go to 10
      if(z.le.2.0d0  ) go to 20
      t=z-2.0d0
      go to 30
10    t=z
      go to 30
20    t=z-1.0d0
30    p=a(18)
      do 40 k1=1,17
      k=18-k1
      p=t*p+a(k)
40    continue
c
      if(z.gt.2.0d0  ) go to 50
      if(z.gt.1.0d0  ) go to 60
      gamfun=p/(z*(z+1.0d0  ))
      return
60    gamfun=p/z
      return
50    gamfun=p
      return
      end
