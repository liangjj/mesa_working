*deck imtql1
      subroutine imtql1 (n, d, e, ierr)
c***begin prologue  imtql1
c***purpose  compute the eigenvalues of a symmetric tridiagonal matrix
c            using the implicit ql method.
c***library   slatec (eispack)
c***category  d4a5, d4c2a
c***type      single precision (imtql1-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure imtql1,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the implicit ql method.
c
c     on input
c
c        n is the order of the matrix.  n is an integer variable.
c
c        d contains the diagonal elements of the symmetric tridiagonal
c          matrix.  d is a one-dimensional real array, dimensioned d(n).
c
c        e contains the subdiagonal elements of the symmetric
c          tridiagonal matrix in its last n-1 positions.  e(1) is
c          arbitrary.  e is a one-dimensional real array, dimensioned
c          e(n).
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an error
c          exit is made, the eigenvalues are correct and ordered for
c          indices 1, 2, ..., ierr-1, but may not be the smallest
c          eigenvalues.
c
c        e has been destroyed.
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c                     the eigenvalues should be correct for indices
c                     1, 2, ..., ierr-1.  these eigenvalues are
c                     ordered, but are not necessarily the smallest.
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
c***end prologue  imtql1
c
      integer i,j,l,m,n,ii,mml,ierr
      real*8 d(*),e(*)
      real*8 b,c,f,g,p,r,s,s1,s2
      real*8 pythag
c
c***first executable statement  imtql1
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      e(n) = 0.0d0
c
      do 290 l = 1, n
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
         if (m .eq. l) go to 215
         if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         g = (d(l+1) - p) / (2.0d0 * e(l))
         r = pythag(g,1.0d0)
         g = d(m) - p + e(l) / (g + sign(r,g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if (abs(f) .lt. abs(g)) go to 150
            c = g / f
            r = sqrt(c*c+1.0d0)
            e(i+1) = f * r
            s = 1.0d0 / r
            c = c * s
            go to 160
  150       s = f / g
            r = sqrt(s*s+1.0d0)
            e(i+1) = g * r
            c = 1.0d0 / r
            s = s * c
  160       g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
  200    continue
c
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
c     .......... order eigenvalues ..........
  215    if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
