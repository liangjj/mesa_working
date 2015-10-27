*deck trbak1
      subroutine trbak1 (nm, n, a, e, m, z)
c***begin prologue  trbak1
c***purpose  form the eigenvectors of real symmetric matrix from
c            the eigenvectors of a symmetric tridiagonal matrix formed
c            by tred1.
c***library   slatec (eispack)
c***category  d4c4
c***type      single precision (trbak1-s)
c***keywords  eigenvectors of a real symmetric matrix, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure trbak1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a real symmetric
c     matrix by back transforming those of the corresponding
c     symmetric tridiagonal matrix determined by  tred1.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, a and z, as declared in the calling
c          program dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix.  n is an integer variable.
c          n must be less than or equal to nm.
c
c        a contains information about the orthogonal transformations
c          used in the reduction by  tred1  in its strict lower
c          triangle.  a is a two-dimensional real array, dimensioned
c          a(nm,n).
c
c        e contains the subdiagonal elements of the tridiagonal matrix
c          in its last n-1 positions.  e(1) is arbitrary.  these
c          elements provide the remaining information about the
c          orthogonal transformations.  e is a one-dimensional real
c          array, dimensioned e(n).
c
c        m is the number of columns of z to be back transformed.
c          m is an integer variable.
c
c        z contains the eigenvectors to be back transformed in its
c          first m columns.  z is a two-dimensional real array,
c          dimensioned z(nm,m).
c
c     on output
c
c        z contains the transformed eigenvectors in its first m columns.
c
c     note that trbak1 preserves vector euclidean norms.
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
c***end prologue  trbak1
c
      integer i,j,k,l,m,n,nm
      real a(nm,*),e(*),z(nm,*)
      real s
c
c***first executable statement  trbak1
      if (m .eq. 0) go to 200
      if (n .eq. 1) go to 200
c
      do 140 i = 2, n
         l = i - 1
         if (e(i) .eq. 0.0e0) go to 140
c
         do 130 j = 1, m
            s = 0.0e0
c
            do 110 k = 1, l
  110       s = s + a(i,k) * z(k,j)
c     .......... divisor below is negative of h formed in tred1.
c                double division avoids possible underflow ..........
            s = (s / a(i,l)) / e(i)
c
            do 120 k = 1, l
  120       z(k,j) = z(k,j) + s * a(i,k)
c
  130    continue
c
  140 continue
c
  200 return
      end
