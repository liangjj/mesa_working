*deck htribk
      subroutine htribk (nm, n, ar, ai, tau, m, zr, zi)
c***begin prologue  htribk
c***purpose  form the eigenvectors of a complex hermitian matrix from
c            the eigenvectors of a real symmetric tridiagonal matrix
c            output from htridi.
c***library   slatec (eispack)
c***category  d4c4
c***type      single precision (htribk-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure trbak1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a complex hermitian
c     matrix by back transforming those of the corresponding
c     real symmetric tridiagonal matrix determined by  htridi.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, ar, ai, zr, and zi, as declared in the
c          calling program dimension statement.  nm is an integer
c          variable.
c
c        n is the order of the matrix.  n is an integer variable.
c          n must be less than or equal to nm.
c
c        ar and ai contain some information about the unitary
c          transformations used in the reduction by  htridi  in the
c          strict lower triangle of ar and the full lower triangle of
c          ai.  the remaining upper parts of the matrices are arbitrary.
c          ar and ai are two-dimensional real arrays, dimensioned
c          ar(nm,n) and ai(nm,n).
c
c        tau contains further information about the transformations.
c          tau is a one-dimensional real array, dimensioned tau(2,n).
c
c        m is the number of eigenvectors to be back transformed.
c          m is an integer variable.
c
c       zr contains the eigenvectors to be back transformed in its first
c          m columns.  the contents of zi are immaterial.  zr and zi are
c          two-dimensional real arrays, dimensioned zr(nm,m) and
c          zi(nm,m).
c
c     on output
c
c        zr and zi contain the real and imaginary parts, respectively,
c          of the transformed eigenvectors in their first m columns.
c
c     note that the last component of each returned vector
c     is real and that vector euclidean norms are preserved.
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
c***end prologue  htribk
c
      integer i,j,k,l,m,n,nm
      real ar(nm,*),ai(nm,*),tau(2,*),zr(nm,*),zi(nm,*)
      real h,s,si
c
c***first executable statement  htribk
      if (m .eq. 0) go to 200
c     .......... transform the eigenvectors of the real symmetric
c                tridiagonal matrix to those of the hermitian
c                tridiagonal matrix. ..........
      do 50 k = 1, n
c
         do 50 j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
   50 continue
c
      if (n .eq. 1) go to 200
c     .......... recover and apply the householder matrices ..........
      do 140 i = 2, n
         l = i - 1
         h = ai(i,i)
         if (h .eq. 0.0e0) go to 140
c
         do 130 j = 1, m
            s = 0.0e0
            si = 0.0e0
c
            do 110 k = 1, l
               s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
               si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
  110       continue
c     .......... double divisions avoid possible underflow ..........
            s = (s / h) / h
            si = (si / h) / h
c
            do 120 k = 1, l
               zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
               zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
  120       continue
c
  130    continue
c
  140 continue
c
  200 return
      end
