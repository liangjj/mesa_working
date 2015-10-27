*deck tred3
      subroutine tred3 (n, nv, a, d, e, e2)
c***begin prologue  tred3
c***purpose  reduce a real symmetric matrix stored in packed form to
c            symmetric tridiagonal matrix using orthogonal
c            transformations.
c***library   slatec (eispack)
c***category  d4c1b1
c***type      single precision (tred3-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure tred3,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix, stored as
c     a one-dimensional array, to a symmetric tridiagonal matrix
c     using orthogonal similarity transformations.
c
c     on input
c
c        n is the order of the matrix a.  n is an integer variable.
c
c        nv is an integer variable set equal to the dimension of the
c          array a as specified in the calling program.  nv must not
c          be less than  n*(n+1)/2.
c
c        a contains the lower triangle, stored row-wise, of the real
c          symmetric packed matrix.  a is a one-dimensional real
c          array, dimensioned a(nv).
c
c     on output
c
c        a contains information about the orthogonal transformations
c          used in the reduction in its first n*(n+1)/2 positions.
c
c        d contains the diagonal elements of the symmetric tridiagonal
c          matrix.  d is a one-dimensional real array, dimensioned d(n).
c
c        e contains the subdiagonal elements of the symmetric
c          tridiagonal matrix in its last n-1 positions.  e(1) is set
c          to zero.  e is a one-dimensional real array, dimensioned
c          e(n).
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c          e2 is a one-dimensional real array, dimensioned e2(n).
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  (none)
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  tred3
c
      integer i,j,k,l,n,ii,iz,jk,nv
      real a(*),d(*),e(*),e2(*)
      real f,g,h,hh,scale
c
c     .......... for i=n step -1 until 1 do -- ..........
c***first executable statement  tred3
      do  300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         iz = (i * l) / 2
         h = 0.0e0
         scale = 0.0e0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
            iz = iz + 1
            d(k) = a(iz)
            scale = scale + abs(d(k))
  120    continue
c
         if (scale .ne. 0.0e0) go to 140
  130    e(i) = 0.0e0
         e2(i) = 0.0e0
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         e2(i) = scale * scale * h
         f = d(l)
         g = -sign(sqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         a(iz) = scale * d(l)
         if (l .eq. 1) go to 290
         f = 0.0e0
c
         do 240 j = 1, l
            g = 0.0e0
            jk = (j * (j-1)) / 2
c     .......... form element of a*u ..........
            do 180 k = 1, l
               jk = jk + 1
               if (k .gt. j) jk = jk + k - 2
               g = g + a(jk) * d(k)
  180       continue
c     .......... form element of p ..........
            e(j) = g / h
            f = f + e(j) * d(j)
  240    continue
c
         hh = f / (h + h)
         jk = 0
c     .......... form reduced a ..........
         do 260 j = 1, l
            f = d(j)
            g = e(j) - hh * f
            e(j) = g
c
            do 260 k = 1, j
               jk = jk + 1
               a(jk) = a(jk) - f * e(k) - g * d(k)
  260    continue
c
  290    d(i) = a(iz+1)
         a(iz+1) = scale * sqrt(h)
  300 continue
c
      return
      end
