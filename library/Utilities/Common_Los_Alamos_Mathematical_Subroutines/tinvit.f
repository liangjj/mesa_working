*deck tinvit
      subroutine tinvit (nm, n, d, e, e2, m, w, ind, z, ierr, rv1, rv2,
     +   rv3, rv4, rv6)
c***begin prologue  tinvit
c***purpose  compute the eigenvectors of symmetric tridiagonal matrix
c            corresponding to specified eigenvalues, using inverse
c            iteration.
c***library   slatec (eispack)
c***category  d4c3
c***type      single precision (tinvit-s)
c***keywords  eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the inverse iteration tech-
c     nique in the algol procedure tristurm by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvectors of a tridiagonal
c     symmetric matrix corresponding to specified eigenvalues,
c     using inverse iteration.
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
c        e2 contains the squares of the corresponding elements of e,
c          with zeros corresponding to negligible elements of e.
c          e(i) is considered negligible if it is not larger than
c          the product of the relative machine precision and the sum
c          of the magnitudes of d(i) and d(i-1).  e2(1) must contain
c          0.0d0 if the eigenvalues are in ascending order, or 2.0d0
c          if the eigenvalues are in descending order.  if  bisect,
c          tridib, or  imtqlv  has been used to find the eigenvalues,
c          their output e2 array is exactly what is expected here.
c          e2 is a one-dimensional real array, dimensioned e2(n).
c
c        m is the number of specified eigenvalues for which eigenvectors
c          are to be determined.  m is an integer variable.
c
c        w contains the m eigenvalues in ascending or descending order.
c          w is a one-dimensional real array, dimensioned w(m).
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc.
c          if  bisect  or  tridib  has been used to determine the
c          eigenvalues, their output ind array is suitable for input
c          to tinvit.  ind is a one-dimensional integer array,
c          dimensioned ind(m).
c
c     on output
c
c       ** all input arrays are unaltered.**
c
c        z contains the associated set of orthonormal eigenvectors.
c          any vector which fails to converge is set to zero.
c          z is a two-dimensional real array, dimensioned z(nm,m).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          -j         if the eigenvector corresponding to the j-th
c                     eigenvalue fails to converge in 5 iterations.
c
c        rv1, rv2 and rv3 are one-dimensional real arrays used for
c          temporary storage.  they are used to store the main diagonal
c          and the two adjacent diagonals of the triangular matrix
c          produced in the inverse iteration process.  rv1, rv2 and
c          rv3 are dimensioned rv1(n), rv2(n) and rv3(n).
c
c        rv4 and rv6 are one-dimensional real arrays used for temporary
c          storage.  rv4 holds the multipliers of the gaussian
c          elimination process.  rv6 holds the approximate eigenvectors
c          in this process.  rv4 and rv6 are dimensioned rv4(n) and
c          rv6(n).
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
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  tinvit
c
      integer i,j,m,n,p,q,r,s,ii,ip,jj,nm,its,tag,ierr,group
      integer ind(*)
      real*8 d(*),e(*),e2(*),w(*),z(nm,*)
      real*8 rv1(*),rv2(*),rv3(*),rv4(*),rv6(*)
      real*8 u,v,uk,xu,x0,x1,eps2,eps3,eps4,norm,order
c
c***first executable statement  tinvit
      ierr = 0
      if (m .eq. 0) go to 1001
      tag = 0
      order = 1.0d0 - e2(1)
      q = 0
c     .......... establish and process next submatrix ..........
  100 p = q + 1
c
      do 120 q = p, n
         if (q .eq. n) go to 140
         if (e2(q+1) .eq. 0.0d0) go to 140
  120 continue
c     .......... find vectors by inverse iteration ..........
  140 tag = tag + 1
      s = 0
c
      do 920 r = 1, m
         if (ind(r) .ne. tag) go to 920
         its = 1
         x1 = w(r)
         if (s .ne. 0) go to 510
c     .......... check for isolated root ..........
         xu = 1.0d0
         if (p .ne. q) go to 490
         rv6(p) = 1.0d0
         go to 870
  490    norm = abs(d(p))
         ip = p + 1
c
         do 500 i = ip, q
  500    norm = max(norm, abs(d(i)) + abs(e(i)))
c     .......... eps2 is the criterion for grouping,
c                eps3 replaces zero pivots and equal
c                roots are modified by eps3,
c                eps4 is taken very small to avoid overflow ..........
         eps2 = 1.0d-3 * norm
         eps3 = norm
  502    eps3 = 0.5d0*eps3
         if (norm + eps3 .gt. norm) go to 502
         uk = sqrt(real(q-p+5))
         eps3 = uk * eps3
         eps4 = uk * eps3
         uk = eps4 / uk
         s = p
  505    group = 0
         go to 520
c     .......... look for close or coincident roots ..........
  510    if (abs(x1-x0) .ge. eps2) go to 505
         group = group + 1
         if (order * (x1 - x0) .le. 0.0d0) x1 = x0 + order * eps3
c     .......... elimination with interchanges and
c                initialization of vector ..........
  520    v = 0.0d0
c
         do 580 i = p, q
            rv6(i) = uk
            if (i .eq. p) go to 560
            if (abs(e(i)) .lt. abs(u)) go to 540
c     .......... warning -- a divide check may occur here if
c                e2 array has not been specified correctly ..........
            xu = u / e(i)
            rv4(i) = xu
            rv1(i-1) = e(i)
            rv2(i-1) = d(i) - x1
            rv3(i-1) = 0.0d0
            if (i .ne. q) rv3(i-1) = e(i+1)
            u = v - xu * rv2(i-1)
            v = -xu * rv3(i-1)
            go to 580
  540       xu = e(i) / u
            rv4(i) = xu
            rv1(i-1) = u
            rv2(i-1) = v
            rv3(i-1) = 0.0d0
  560       u = d(i) - x1 - xu * v
            if (i .ne. q) v = e(i+1)
  580    continue
c
         if (u .eq. 0.0d0) u = eps3
         rv1(q) = u
         rv2(q) = 0.0d0
         rv3(q) = 0.0d0
c     .......... back substitution
c                for i=q step -1 until p do -- ..........
  600    do 620 ii = p, q
            i = p + q - ii
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
            v = u
            u = rv6(i)
  620    continue
c     .......... orthogonalize with respect to previous
c                members of group ..........
         if (group .eq. 0) go to 700
         j = r
c
         do 680 jj = 1, group
  630       j = j - 1
            if (ind(j) .ne. tag) go to 630
            xu = 0.0d0
c
            do 640 i = p, q
  640       xu = xu + rv6(i) * z(i,j)
c
            do 660 i = p, q
  660       rv6(i) = rv6(i) - xu * z(i,j)
c
  680    continue
c
  700    norm = 0.0d0
c
         do 720 i = p, q
  720    norm = norm + abs(rv6(i))
c
         if (norm .ge. 1.0d0) go to 840
c     .......... forward substitution ..........
         if (its .eq. 5) go to 830
         if (norm .ne. 0.0d0) go to 740
         rv6(s) = eps4
         s = s + 1
         if (s .gt. q) s = p
         go to 780
  740    xu = eps4 / norm
c
         do 760 i = p, q
  760    rv6(i) = rv6(i) * xu
c     .......... elimination operations on next vector
c                iterate ..........
  780    do 820 i = ip, q
            u = rv6(i)
c     .......... if rv1(i-1) .eq. e(i), a row interchange
c                was performed earlier in the
c                triangularization process ..........
            if (rv1(i-1) .ne. e(i)) go to 800
            u = rv6(i-1)
            rv6(i-1) = rv6(i)
  800       rv6(i) = u - rv4(i) * rv6(i-1)
  820    continue
c
         its = its + 1
         go to 600
c     .......... set error -- non-converged eigenvector ..........
  830    ierr = -r
         xu = 0.0d0
         go to 870
c     .......... normalize so that sum of squares is
c                1 and expand to full order ..........
  840    u = 0.0d0
c
         do 860 i = p, q
  860    u = u + rv6(i)**2
c
         xu = 1.0d0 / sqrt(u)
c
  870    do 880 i = 1, n
  880    z(i,r) = 0.0d0
c
         do 900 i = p, q
  900    z(i,r) = rv6(i) * xu
c
         x0 = x1
  920 continue
c
      if (q .lt. n) go to 100
 1001 return
      end
