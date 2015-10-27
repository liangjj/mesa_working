*deck trbak3
      subroutine trbak3 (nm, n, nv, a, m, z)
c***begin prologue  trbak3
c***purpose  form the eigenvectors of a real symmetric matrix from the
c            eigenvectors of a symmetric tridiagonal matrix formed
c            by tred3.
c***library   slatec (eispack)
c***category  d4c4
c***type      single precision (trbak3-s)
c***keywords  eigenvectors of a real symmetric matrix, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of the algol procedure trbak3,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a real symmetric
c     matrix by back transforming those of the corresponding
c     symmetric tridiagonal matrix determined by  tred3.
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
c        nv is an integer variable set equal to the dimension of the
c          array a as specified in the calling program.  nv must not
c          be less than  n*(n+1)/2.
c
c        a contains information about the orthogonal transformations
c          used in the reduction by  tred3  in its first n*(n+1)/2
c          positions.  a is a one-dimensional real array, dimensioned
c          a(nv).
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
c     note that trbak3 preserves vector euclidean norms.
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
c***end prologue  trbak3
c
      integer i,j,k,l,m,n,ik,iz,nm,nv
      real a(*),z(nm,*)
      real h,s
c
c***first executable statement  trbak3
      if (m .eq. 0) go to 200
      if (n .eq. 1) go to 200
c
      do 140 i = 2, n
         l = i - 1
         iz = (i * l) / 2
         ik = iz + i
         h = a(ik)
         if (h .eq. 0.0e0) go to 140
c
         do 130 j = 1, m
            s = 0.0e0
            ik = iz
c
            do 110 k = 1, l
               ik = ik + 1
               s = s + a(ik) * z(k,j)
  110       continue
c     .......... double division avoids possible underflow ..........
            s = (s / h) / h
            ik = iz
c
            do 120 k = 1, l
               ik = ik + 1
               z(k,j) = z(k,j) - s * a(ik)
  120       continue
c
  130    continue
c
  140 continue
c
  200 return
      end
