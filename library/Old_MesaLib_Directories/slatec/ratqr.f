*deck ratqr
      subroutine ratqr (n, eps1, d, e, e2, m, w, ind, bd, type, idef,
     +   ierr)
c***begin prologue  ratqr
c***purpose  compute the largest or smallest eigenvalues of a symmetric
c            tridiagonal matrix using the rational qr method with newton
c            correction.
c***library   slatec (eispack)
c***category  d4a5, d4c2a
c***type      single precision (ratqr-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure ratqr,
c     num. math. 11, 264-272(1968) by reinsch and bauer.
c     handbook for auto. comp., vol.ii-linear algebra, 257-265(1971).
c
c     this subroutine finds the algebraically smallest or largest
c     eigenvalues of a symmetric tridiagonal matrix by the
c     rational qr method with newton corrections.
c
c     on input
c
c        n is the order of the matrix.  n is an integer variable.
c
c        eps1 is a theoretical absolute error tolerance for the
c          computed eigenvalues.  if the input eps1 is non-positive, or
c          indeed smaller than its default value, it is reset at each
c          iteration to the respective default value, namely, the
c          product of the relative machine precision and the magnitude
c          of the current eigenvalue iterate.  the theoretical absolute
c          error in the k-th eigenvalue is usually not greater than
c          k times eps1.  eps1 is a real variable.
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
c        m is the number of eigenvalues to be found.  m is an integer
c          variable.
c
c        idef should be set to 1 if the input matrix is known to be
c          positive definite, to -1 if the input matrix is known to
c          be negative definite, and to 0 otherwise.  idef is an
c          integer variable.
c
c        type should be set to .true. if the smallest eigenvalues are
c          to be found, and to .false. if the largest eigenvalues are
c          to be found.  type is a logical variable.
c
c     on output
c
c        eps1 is unaltered unless it has been reset to its
c          (last) default value.
c
c        d and e are unaltered (unless w overwrites d).
c
c        elements of e2, corresponding to elements of e regarded as
c          negligible, have been replaced by zero causing the matrix
c          to split into a direct sum of submatrices.  e2(1) is set
c          to 0.0e0 if the smallest eigenvalues have been found, and
c          to 2.0e0 if the largest eigenvalues have been found.  e2
c          is otherwise unaltered (unless overwritten by bd).
c
c        w contains the m algebraically smallest eigenvalues in
c          ascending order, or the m largest eigenvalues in descending
c          order.  if an error exit is made because of an incorrect
c          specification of idef, no eigenvalues are found.  if the
c          newton iterates for a particular eigenvalue are not monotone,
c          the best estimate obtained is returned and ierr is set.
c          w is a one-dimensional real array, dimensioned w(n).  w need
c          not be distinct from d.
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc.
c          ind is an one-dimensional integer array, dimensioned ind(n).
c
c        bd contains refined bounds for the theoretical errors of the
c          corresponding eigenvalues in w.  these bounds are usually
c          within the tolerance specified by eps1.  bd is a one-
c          dimensional real array, dimensioned bd(n).  bd need not be
c          distinct from e2.
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          6*n+1      if  idef  is set to 1 and  type  to .true.
c                     when the matrix is not positive definite, or
c                     if  idef  is set to -1 and  type  to .false.
c                     when the matrix is not negative definite,
c                     no eigenvalues are computed, or
c                     m is greater than n,
c          5*n+k      if successive iterates to the k-th eigenvalue
c                     are not monotone increasing, where k refers
c                     to the last such occurrence.
c
c     note that subroutine tridib is generally faster and more
c     accurate than ratqr if the eigenvalues are clustered.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  r1mach
c***revision history  (yymmdd)
c   760101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  ratqr
c
      integer i,j,k,m,n,ii,jj,k1,idef,ierr,jdef
      real d(*),e(*),e2(*),w(*),bd(*)
      real f,p,q,r,s,ep,qp,err,tot,eps1,delta,machep
      integer ind(*)
      logical first, type
c
      save first, machep
      data first /.true./
c***first executable statement  ratqr
      if (first) then
         machep = r1mach(4)
      endif
      first = .false.
c
      ierr = 0
      jdef = idef
c     .......... copy d array into w ..........
      do 20 i = 1, n
   20 w(i) = d(i)
c
      if (type) go to 40
      j = 1
      go to 400
   40 err = 0.0e0
      s = 0.0e0
c     .......... look for small sub-diagonal entries and define
c                initial shift from lower gerschgorin bound.
c                copy e2 array into bd ..........
      tot = w(1)
      q = 0.0e0
      j = 0
c
      do 100 i = 1, n
         p = q
         if (i .eq. 1) go to 60
         if (p .gt. machep * (abs(d(i)) + abs(d(i-1)))) go to 80
   60    e2(i) = 0.0e0
   80    bd(i) = e2(i)
c     .......... count also if element of e2 has underflowed ..........
         if (e2(i) .eq. 0.0e0) j = j + 1
         ind(i) = j
         q = 0.0e0
         if (i .ne. n) q = abs(e(i+1))
         tot = min(w(i)-p-q,tot)
  100 continue
c
      if (jdef .eq. 1 .and. tot .lt. 0.0e0) go to 140
c
      do 110 i = 1, n
  110 w(i) = w(i) - tot
c
      go to 160
  140 tot = 0.0e0
c
  160 do 360 k = 1, m
c     .......... next qr transformation ..........
  180    tot = tot + s
         delta = w(n) - s
         i = n
         f = abs(machep*tot)
         if (eps1 .lt. f) eps1 = f
         if (delta .gt. eps1) go to 190
         if (delta .lt. (-eps1)) go to 1000
         go to 300
c     .......... replace small sub-diagonal squares by zero
c                to reduce the incidence of underflows ..........
  190    if (k .eq. n) go to 210
         k1 = k + 1
         do 200 j = k1, n
            if (bd(j) .le. (machep*(w(j)+w(j-1))) ** 2) bd(j) = 0.0e0
  200    continue
c
  210    f = bd(n) / delta
         qp = delta + f
         p = 1.0e0
         if (k .eq. n) go to 260
         k1 = n - k
c     .......... for i=n-1 step -1 until k do -- ..........
         do 240 ii = 1, k1
            i = n - ii
            q = w(i) - s - f
            r = q / qp
            p = p * r + 1.0e0
            ep = f * r
            w(i+1) = qp + ep
            delta = q - ep
            if (delta .gt. eps1) go to 220
            if (delta .lt. (-eps1)) go to 1000
            go to 300
  220       f = bd(i) / q
            qp = delta + f
            bd(i+1) = qp * ep
  240    continue
c
  260    w(k) = qp
         s = qp / p
         if (tot + s .gt. tot) go to 180
c     .......... set error -- irregular end of iteration.
c                deflate minimum diagonal element ..........
         ierr = 5 * n + k
         s = 0.0e0
         delta = qp
c
         do 280 j = k, n
            if (w(j) .gt. delta) go to 280
            i = j
            delta = w(j)
  280    continue
c     .......... convergence ..........
  300    if (i .lt. n) bd(i+1) = bd(i) * f / qp
         ii = ind(i)
         if (i .eq. k) go to 340
         k1 = i - k
c     .......... for j=i-1 step -1 until k do -- ..........
         do 320 jj = 1, k1
            j = i - jj
            w(j+1) = w(j) - s
            bd(j+1) = bd(j)
            ind(j+1) = ind(j)
  320    continue
c
  340    w(k) = tot
         err = err + abs(delta)
         bd(k) = err
         ind(k) = ii
  360 continue
c
      if (type) go to 1001
      f = bd(1)
      e2(1) = 2.0e0
      bd(1) = f
      j = 2
c     .......... negate elements of w for largest values ..........
  400 do 500 i = 1, n
  500 w(i) = -w(i)
c
      jdef = -jdef
      go to (40,1001), j
c     .......... set error -- idef specified incorrectly ..........
 1000 ierr = 6 * n + 1
 1001 return
      end
