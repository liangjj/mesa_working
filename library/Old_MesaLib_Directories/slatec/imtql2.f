*deck imtql2
      subroutine imtql2 (nm, n, d, e, z, ierr)
c***begin prologue  imtql2
c***purpose  compute the eigenvalues and eigenvectors of a symmetric
c            tridiagonal matrix using the implicit ql method.
c***library   slatec (eispack)
c***category  d4a5, d4c2a
c***type      single precision (imtql2-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure imtql2,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the implicit ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameter, z, as declared in the calling program
c          dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix.  n is an integer variable.
c          n must be less than or equal to nm.
c
c        d contains the diagonal elements of the symmetric tridiagonal
c          matrix.  d is a one-dimensional real array, dimensioned d(n).
c
c        e contains the subdiagonal elements of the symmetric
c          tridiagonal matrix in its last n-1 positions.  e(1) is
c          arbitrary.  e is a one-dimensional real array, dimensioned
c          e(n).
c
c        z contains the transformation matrix produced in the reduction
c          by  tred2,  if performed.  this transformation matrix is
c          necessary if you want to obtain the eigenvectors of the full
c          symmetric matrix.  if the eigenvectors of the symmetric
c          tridiagonal matrix are desired, z must contain the identity
c          matrix.  z is a two-dimensional real array, dimensioned
c          z(nm,n).
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1, 2, ..., ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the full symmetric
c          or symmetric tridiagonal matrix, depending on what it
c          contained on input.  if an error exit is made,  z contains
c          the eigenvectors associated with the stored eigenvalues.
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c                     the eigenvalues and eigenvectors should be correct
c                     for indices 1, 2, ..., ierr-1, but the eigenvalues
c                     are not ordered.
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
c***end prologue  imtql2
c
      integer i,j,k,l,m,n,ii,nm,mml,ierr
      real d(*),e(*),z(nm,*)
      real b,c,f,g,p,r,s,s1,s2
      real pythag
c
c***first executable statement  imtql2
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      e(n) = 0.0e0
c
      do 240 l = 1, n
         j = 0
c     .......... look for small sub-diagonal element ..........
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            s1 = abs(d(m)) + abs(d(m+1))
            s2 = s1 + abs(e(m))
            if (s2 .eq. s1) go to 120
  110    continue
c
  120    p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         g = (d(l+1) - p) / (2.0e0 * e(l))
         r = pythag(g,1.0e0)
         g = d(m) - p + e(l) / (g + sign(r,g))
         s = 1.0e0
         c = 1.0e0
         p = 0.0e0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if (abs(f) .lt. abs(g)) go to 150
            c = g / f
            r = sqrt(c*c+1.0e0)
            e(i+1) = f * r
            s = 1.0e0 / r
            c = c * s
            go to 160
  150       s = f / g
            r = sqrt(s*s+1.0e0)
            e(i+1) = g * r
            c = 1.0e0 / r
            s = s * c
  160       g = d(i+1) - p
            r = (d(i) - g) * s + 2.0e0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
c     .......... form vector ..........
            do 180 k = 1, n
               f = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * f
               z(k,i) = c * z(k,i) - s * f
  180       continue
c
  200    continue
c
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0e0
         go to 105
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
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
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
