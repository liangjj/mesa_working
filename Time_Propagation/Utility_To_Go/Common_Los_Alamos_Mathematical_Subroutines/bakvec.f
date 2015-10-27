*deck bakvec
      subroutine bakvec (nm, n, t, e, m, z, ierr)
c***begin prologue  bakvec
c***purpose  form the eigenvectors of a certain real non-symmetric
c            tridiagonal matrix from a symmetric tridiagonal matrix
c            output from figi.
c***library   slatec (eispack)
c***category  d4c4
c***type      single precision (bakvec-s)
c***keywords  eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine forms the eigenvectors of a nonsymmetric
c     tridiagonal matrix by back transforming those of the
c     corresponding symmetric matrix determined by  figi.
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
c        e contains the subdiagonal elements of the symmetric
c          matrix in its last n-1 positions.  e(1) is arbitrary.
c          e is a one-dimensional real array, dimensioned e(n).
c
c        m is the number of eigenvectors to be back transformed.
c          m is an integer variable.
c
c        z contains the eigenvectors to be back transformed
c          in its first m columns.  z is a two-dimensional real
c          array, dimensioned z(nm,m).
c
c     on output
c
c        t is unaltered.
c
c        e is destroyed.
c
c        z contains the transformed eigenvectors in its first m columns.
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          2*n+i      if e(i) is zero with t(i,1) or t(i-1,3) non-zero.
c                     in this case, the symmetric matrix is not similar
c                     to the original matrix, and the eigenvectors
c                     cannot be found by this program.
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
c***end prologue  bakvec
c
      integer i,j,m,n,nm,ierr
      real*8 t(nm,3),e(*),z(nm,*)
c
c***first executable statement  bakvec
      ierr = 0
      if (m .eq. 0) go to 1001
      e(1) = 1.0d0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
         if (e(i) .ne. 0.0d0) go to 80
         if (t(i,1) .ne. 0.0d0 .or. t(i-1,3) .ne. 0.0d0) go to 1000
         e(i) = 1.0d0
         go to 100
   80    e(i) = e(i-1) * e(i) / t(i-1,3)
  100 continue
c
      do 120 j = 1, m
c
         do 120 i = 2, n
         z(i,j) = z(i,j) * e(i)
  120 continue
c
      go to 1001
c     .......... set error -- eigenvectors cannot be
c                found by this program ..........
 1000 ierr = 2 * n + i
 1001 return
      end
