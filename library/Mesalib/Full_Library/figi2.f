*deck figi2
      subroutine figi2 (nm, n, t, d, e, z, ierr)
c***begin prologue  figi2
c***purpose  transforms certain real non-symmetric tridiagonal matrix
c            to symmetric tridiagonal matrix.
c***library   slatec (eispack)
c***category  d4c1c
c***type      single precision (figi2-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     given a nonsymmetric tridiagonal matrix such that the products
c     of corresponding pairs of off-diagonal elements are all
c     non-negative, and zero only when both factors are zero, this
c     subroutine reduces it to a symmetric tridiagonal matrix
c     using and accumulating diagonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, t and z, as declared in the calling
c          program dimension statement.  nm is an integer variable.
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
c        z contains the diagonal transformation matrix produced in the
c          symmetrization.  z is a two-dimensional real array,
c          dimensioned z(nm,n).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          n+i        if t(i,1)*t(i-1,3) is negative,
c          2*n+i      if t(i,1)*t(i-1,3) is zero with one factor
c                     non-zero.  in these cases, there does not exist
c                     a symmetrizing similarity transformation which
c                     is essential for the validity of the later
c                     eigenvector computation.
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
c***end prologue  figi2
c
      integer i,j,n,nm,ierr
      real*8 t(nm,3),d(*),e(*),z(nm,*)
      real*8 h
c
c***first executable statement  figi2
      ierr = 0
c
      do 100 i = 1, n
c
         do 50 j = 1, n
   50    z(i,j) = 0.0d0
c
         if (i .eq. 1) go to 70
         h = t(i,1) * t(i-1,3)
         if (h) 900, 60, 80
   60    if (t(i,1) .ne. 0.0d0 .or. t(i-1,3) .ne. 0.0d0) go to 1000
         e(i) = 0.0d0
   70    z(i,i) = 1.0d0
         go to 90
   80    e(i) = sqrt(h)
         z(i,i) = z(i-1,i-1) * e(i) / t(i-1,3)
   90    d(i) = t(i,2)
  100 continue
c
      go to 1001
c     .......... set error -- product of some pair of off-diagonal
c                elements is negative ..........
  900 ierr = n + i
      go to 1001
c     .......... set error -- product of some pair of off-diagonal
c                elements is zero with one member non-zero ..........
 1000 ierr = 2 * n + i
 1001 return
      end
