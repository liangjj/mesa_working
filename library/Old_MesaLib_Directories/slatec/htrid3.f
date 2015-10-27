*deck htrid3
      subroutine htrid3 (nm, n, a, d, e, e2, tau)
c***begin prologue  htrid3
c***purpose  reduce a complex hermitian (packed) matrix to a real
c            symmetric tridiagonal matrix by unitary similarity
c            transformations.
c***library   slatec (eispack)
c***category  d4c1b1
c***type      single precision (htrid3-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure tred3, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a complex hermitian matrix, stored as
c     a single square array, to a real symmetric tridiagonal matrix
c     using unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameter, a, as declared in the calling program
c          dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix.  n is an integer variable.
c          n must be less than or equal to nm.
c
c        a contains the lower triangle of the complex hermitian input
c          matrix.  the real parts of the matrix elements are stored
c          in the full lower triangle of a, and the imaginary parts
c          are stored in the transposed positions of the strict upper
c          triangle of a.  no storage is required for the zero
c          imaginary parts of the diagonal elements.  a is a two-
c          dimensional real array, dimensioned a(nm,n).
c
c     on output
c
c        a contains some information about the unitary transformations
c          used in the reduction.
c
c        d contains the diagonal elements of the real symmetric
c          tridiagonal matrix.  d is a one-dimensional real array,
c          dimensioned d(n).
c
c        e contains the subdiagonal elements of the real tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c          e is a one-dimensional real array, dimensioned e(n).
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is set to zero.  e2 may coincide with e if the squares
c          are not needed.  e2 is a one-dimensional real array,
c          dimensioned e2(n).
c
c        tau contains further information about the transformations.
c          tau is a one-dimensional real array, dimensioned tau(2,n).
c
c     calls pythag(a,b) for sqrt(a**2 + b**2).
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  pythag
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  htrid3
c
      integer i,j,k,l,n,ii,nm,jm1,jp1
      real a(nm,*),d(*),e(*),e2(*),tau(2,*)
      real f,g,h,fi,gi,hh,si,scale
      real pythag
c
c***first executable statement  htrid3
      tau(1,n) = 1.0e0
      tau(2,n) = 0.0e0
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0e0
         scale = 0.0e0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + abs(a(i,k)) + abs(a(k,i))
c
         if (scale .ne. 0.0e0) go to 140
         tau(1,l) = 1.0e0
         tau(2,l) = 0.0e0
  130    e(i) = 0.0e0
         e2(i) = 0.0e0
         go to 290
c
  140    do 150 k = 1, l
            a(i,k) = a(i,k) / scale
            a(k,i) = a(k,i) / scale
            h = h + a(i,k) * a(i,k) + a(k,i) * a(k,i)
  150    continue
c
         e2(i) = scale * scale * h
         g = sqrt(h)
         e(i) = scale * g
         f = pythag(a(i,l),a(l,i))
c     .......... form next diagonal element of matrix t ..........
         if (f .eq. 0.0e0) go to 160
         tau(1,l) = (a(l,i) * tau(2,i) - a(i,l) * tau(1,i)) / f
         si = (a(i,l) * tau(2,i) + a(l,i) * tau(1,i)) / f
         h = h + f * g
         g = 1.0e0 + g / f
         a(i,l) = g * a(i,l)
         a(l,i) = g * a(l,i)
         if (l .eq. 1) go to 270
         go to 170
  160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         a(i,l) = g
  170    f = 0.0e0
c
         do 240 j = 1, l
            g = 0.0e0
            gi = 0.0e0
            if (j .eq. 1) go to 190
            jm1 = j - 1
c     .......... form element of a*u ..........
            do 180 k = 1, jm1
               g = g + a(j,k) * a(i,k) + a(k,j) * a(k,i)
               gi = gi - a(j,k) * a(k,i) + a(k,j) * a(i,k)
  180       continue
c
  190       g = g + a(j,j) * a(i,j)
            gi = gi - a(j,j) * a(j,i)
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + a(k,j) * a(i,k) - a(j,k) * a(k,i)
               gi = gi - a(k,j) * a(k,i) - a(j,k) * a(i,k)
  200       continue
c     .......... form element of p ..........
  220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * a(i,j) - tau(2,j) * a(j,i)
  240    continue
c
         hh = f / (h + h)
c     .......... form reduced a ..........
         do 260 j = 1, l
            f = a(i,j)
            g = e(j) - hh * f
            e(j) = g
            fi = -a(j,i)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
            a(j,j) = a(j,j) - 2.0e0 * (f * g + fi * gi)
            if (j .eq. 1) go to 260
            jm1 = j - 1
c
            do 250 k = 1, jm1
               a(j,k) = a(j,k) - f * e(k) - g * a(i,k)
     1                         + fi * tau(2,k) + gi * a(k,i)
               a(k,j) = a(k,j) - f * tau(2,k) - g * a(k,i)
     1                         - fi * e(k) - gi * a(i,k)
  250       continue
c
  260    continue
c
  270    do 280 k = 1, l
            a(i,k) = scale * a(i,k)
            a(k,i) = scale * a(k,i)
  280    continue
c
         tau(2,l) = -si
  290    d(i) = a(i,i)
         a(i,i) = scale * sqrt(h)
  300 continue
c
      return
      end
