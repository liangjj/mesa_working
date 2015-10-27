*deck figi
      subroutine figi (nm, n, t, d, e, e2, ierr)
c***begin prologue  figi
c***purpose  transforms certain real non-symmetric tridiagonal matrix
c            to symmetric tridiagonal matrix.
c***library   slatec (eispack)
c***category  d4c1c
c***type      single precision (figi-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     given a nonsymmetric tridiagonal matrix such that the products
c     of corresponding pairs of off-diagonal elements are all
c     non-negative, this subroutine reduces it to a symmetric
c     tridiagonal matrix with the same eigenvalues.  if, further,
c     a zero product only occurs when both factors are zero,
c     the reduced matrix is similar to the original matrix.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameter, t, as declared in the calling program
c          dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix t.  n is an integer variable.
c          n must be less than or equal to nm.
c
c        t contains the nonsymmetric matrix.  its subdiagonal is
c          stored in the last n-1 positions of the first column,
c          its diagonal in the n positions of the second column,
c          and its superdiagonal in the first n-1 positions of
c          the third column.  t(1,1) and t(n,3) are arbitrary.
c          t is a two-dimensional real array, dimensioned t(nm,3).
c
c     on output
c
c        t is unaltered.
c
c        d contains the diagonal elements of the tridiagonal symmetric
c          matrix.  d is a one-dimensional real array, dimensioned d(n).
c
c        e contains the subdiagonal elements of the tridiagonal
c          symmetric matrix in its last n-1 positions.  e(1) is not set.
c          e is a one-dimensional real array, dimensioned e(n).
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c          e2 is a one-dimensional real array, dimensioned e2(n).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          n+i        if t(i,1)*t(i-1,3) is negative and a symmetric
c                     matrix cannot be produced with figi,
c          -(3*n+i)   if t(i,1)*t(i-1,3) is zero with one factor
c                     non-zero.  in this case, the eigenvectors of
c                     the symmetric matrix are not simply related
c                     to those of  t  and should not be sought.
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
c***end prologue  figi
c
      integer i,n,nm,ierr
      real t(nm,3),d(*),e(*),e2(*)
c
c***first executable statement  figi
      ierr = 0
c
      do 100 i = 1, n
         if (i .eq. 1) go to 90
         e2(i) = t(i,1) * t(i-1,3)
         if (e2(i)) 1000, 60, 80
   60    if (t(i,1) .eq. 0.0e0 .and. t(i-1,3) .eq. 0.0e0) go to 80
c     .......... set error -- product of some pair of off-diagonal
c                elements is zero with one member non-zero ..........
         ierr = -(3 * n + i)
   80    e(i) = sqrt(e2(i))
   90    d(i) = t(i,2)
  100 continue
c
      go to 1001
c     .......... set error -- product of some pair of off-diagonal
c                elements is negative ..........
 1000 ierr = n + i
 1001 return
      end
