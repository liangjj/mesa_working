*deck imtqlv
      subroutine imtqlv (n, d, e, e2, w, ind, ierr, rv1)
c***begin prologue  imtqlv
c***purpose  compute the eigenvalues of a symmetric tridiagonal matrix
c            using the implicit ql method.  eigenvectors may be computed
c            later.
c***library   slatec (eispack)
c***category  d4a5, d4c2a
c***type      single precision (imtqlv-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a variant of  imtql1  which is a translation of
c     algol procedure imtql1, num. math. 12, 377-383(1968) by martin and
c     wilkinson, as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues of a symmetric tridiagonal
c     matrix by the implicit ql method and associates with them
c     their corresponding submatrix indices.
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
c        e2 contains the squares of the corresponding elements of e in
c          its last n-1 positions.  e2(1) is arbitrary.  e2 is a one-
c          dimensional real array, dimensioned e2(n).
c
c     on output
c
c        d and e are unaltered.
c
c        elements of e2, corresponding to elements of e regarded as
c          negligible, have been replaced by zero causing the matrix to
c          split into a direct sum of submatrices.  e2(1) is also set
c          to zero.
c
c        w contains the eigenvalues in ascending order.  if an error
c          exit is made, the eigenvalues are correct and ordered for
c          indices 1, 2, ..., ierr-1, but may not be the smallest
c          eigenvalues.  w is a one-dimensional real array, dimensioned
c          w(n).
c
c        ind contains the submatrix indices associated with the
c          corresponding eigenvalues in w -- 1 for eigenvalues belonging
c          to the first submatrix from the top, 2 for those belonging to
c          the second submatrix, etc.  ind is a one-dimensional real
c          array, dimensioned ind(n).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c                     the eigenvalues should be correct for indices
c                     1, 2, ..., ierr-1.  these eigenvalues are
c                     ordered, but are not necessarily the smallest.
c
c        rv1 is a one-dimensional real array used for temporary storage,
c          dimensioned rv1(n).
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
c***end prologue  imtqlv
c
      integer i,j,k,l,m,n,ii,mml,tag,ierr
      real d(*),e(*),e2(*),w(*),rv1(*)
      real b,c,f,g,p,r,s,s1,s2
      real pythag
      integer ind(*)
c
c***first executable statement  imtqlv
      ierr = 0
      k = 0
      tag = 0
c
      do 100 i = 1, n
         w(i) = d(i)
         if (i .ne. 1) rv1(i-1) = e(i)
  100 continue
c
      e2(1) = 0.0e0
      rv1(n) = 0.0e0
c
      do 290 l = 1, n
         j = 0
c     .......... look for small sub-diagonal element ..........
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            s1 = abs(w(m)) + abs(w(m+1))
            s2 = s1 + abs(rv1(m))
            if (s2 .eq. s1) go to 120
c     .......... guard against underflowed element of e2 ..........
            if (e2(m+1) .eq. 0.0e0) go to 125
  110    continue
c
  120    if (m .le. k) go to 130
         if (m .ne. n) e2(m+1) = 0.0e0
  125    k = m
         tag = tag + 1
  130    p = w(l)
         if (m .eq. l) go to 215
         if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         g = (w(l+1) - p) / (2.0e0 * rv1(l))
         r = pythag(g,1.0e0)
         g = w(m) - p + rv1(l) / (g + sign(r,g))
         s = 1.0e0
         c = 1.0e0
         p = 0.0e0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            f = s * rv1(i)
            b = c * rv1(i)
            if (abs(f) .lt. abs(g)) go to 150
            c = g / f
            r = sqrt(c*c+1.0e0)
            rv1(i+1) = f * r
            s = 1.0e0 / r
            c = c * s
            go to 160
  150       s = f / g
            r = sqrt(s*s+1.0e0)
            rv1(i+1) = g * r
            c = 1.0e0 / r
            s = s * c
  160       g = w(i+1) - p
            r = (w(i) - g) * s + 2.0e0 * c * b
            p = s * r
            w(i+1) = g + p
            g = c * r - b
  200    continue
c
         w(l) = w(l) - p
         rv1(l) = g
         rv1(m) = 0.0e0
         go to 105
c     .......... order eigenvalues ..........
  215    if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. w(i-1)) go to 270
            w(i) = w(i-1)
            ind(i) = ind(i-1)
  230    continue
c
  250    i = 1
  270    w(i) = p
         ind(i) = tag
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
