*deck eisdoc
      subroutine eisdoc
c***begin prologue  eisdoc
c***purpose  documentation for eispack, a collection of subprograms for
c            solving matrix eigen-problems.
c***library   slatec (eispack)
c***category  d4, z
c***type      all (eisdoc-a)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  vandevender, w. h., (snla)
c***description
c
c                 **********eispack routines**********
c
c single double complx
c ------ ------ ------
c
c rs       -    ch     computes eigenvalues and, optionally,
c                      eigenvectors of real symmetric
c                      (complex hermitian) matrix.
c
c rsp      -      -    compute eigenvalues and, optionally,
c                      eigenvectors of real symmetric matrix
c                      packed into a one dimensional array.
c
c rg       -    cg     computes eigenvalues and, optionally,
c                      eigenvectors of a real (complex) general
c                      matrix.
c
c bisect   -      -    compute eigenvalues of symmetric tridiagonal
c                      matrix given interval using sturm sequencing.
c
c imtql1   -      -    computes eigenvalues of symmetric tridiagonal
c                      matrix implicit ql method.
c
c imtql2   -      -    computes eigenvalues and eigenvectors of
c                      symmetric tridiagonal matrix using
c                      implicit ql method.
c
c imtqlv   -      -    computes eigenvalues of symmetric tridiagonal
c                      matrix by the implicit ql method.
c                      eigenvectors may be computed later.
c
c ratqr    -      -    computes largest or smallest eigenvalues
c                      of symmetric tridiagonal matrix using
c                      rational qr method with newton correction.
c
c rst      -      -    compute eigenvalues and, optionally,
c                      eigenvectors of real symmetric tridiagonal
c                      matrix.
c
c rt       -      -    compute eigenvalues and eigenvectors of
c                      a special real tridiagonal matrix.
c
c tql1     -      -    compute eigenvalues of symmetric tridiagonal
c                      matrix by ql method.
c
c tql2     -      -    compute eigenvalues and eigenvectors
c                      of symmetric tridiagonal matrix.
c
c tqlrat   -      -    computes eigenvalues of symmetric
c                      tridiagonal matrix a rational variant
c                      of the ql method.
c
c tridib   -      -    computes eigenvalues of symmetric
c                      tridiagonal matrix given interval using
c                      sturm sequencing.
c
c tsturm   -      -    computes eigenvalues of symmetric tridiagonal
c                      matrix given interval and eigenvectors
c                      by sturm sequencing.  this subroutine
c                      is a translation of the algol procedure
c                      tristurm by peters and wilkinson. handbook
c                      for auto. comp., vol.ii-linear algebra,
c                      418-439(1971).
c
c bqr      -      -    computes some of the eigenvalues of a real
c                      symmetric matrix using the qr method with
c                      shifts of origin.
c
c rsb      -      -    computes eigenvalues and, optionally,
c                      eigenvectors of symmetric band matrix.
c
c rsg      -      -    computes eigenvalues and, optionally,
c                      eigenvectors of symmetric generalized
c                      eigenproblem: a*x=(lambda)*b*x
c
c rsgab    -      -    computes eigenvalues and, optionally,
c                      eigenvectors of symmetric generalized
c                      eigenproblem: a*b*x=(lambda)*x
c
c rsgba    -      -    computes eigenvalues and, optionally,
c                      eigenvectors of symmetric generalized
c                      eigenproblem: b*a*x=(lambda)*x
c
c rgg      -      -    computes eigenvalues and eigenvectors
c                      for real generalized eigenproblem:
c                      a*x=(lambda)*b*x.
c
c balanc   -    cbal   balances a general real (complex)
c                      matrix and isolates eigenvalues whenever
c                      possible.
c
c bandr    -      -    reduces real symmetric band matrix
c                      to symmetric tridiagonal matrix and,
c                      optionally, accumulates orthogonal similarity
c                      transformations.
c
c htrid3   -      -    reduces complex hermitian (packed) matrix
c                      to real symmetric tridiagonal matrix by unitary
c                      similarity transformations.
c
c htridi   -      -    reduces complex hermitian matrix to real
c                      symmetric tridiagonal matrix using unitary
c                      similarity transformations.
c
c tred1    -      -    reduce real symmetric matrix to symmetric
c                      tridiagonal matrix using orthogonal
c                      similarity transformations.
c
c tred2    -      -    reduce real symmetric matrix to symmetric
c                      tridiagonal matrix using and accumulating
c                      orthogonal transformations.
c
c tred3    -      -    reduce  symmetric matrix stored in packed
c                      form to symmetric tridiagonal matrix using
c                      orthogonal transformations.
c
c elmhes   -    comhes reduces real (complex) general matrix to
c                      upper hessenberg form using stabilized
c                      elementary similarity transformations.
c
c orthes   -    corth  reduces real (complex) general matrix to upper
c                      hessenberg form orthogonal (unitary)
c                      similarity transformations.
c
c qzhes    -      -    the first step of the qz algorithm for solving
c                      generalized matrix eigenproblems.  accepts
c                      a pair of real general matrices and reduces
c                      one of them to upper hessenberg and the other
c                      to upper triangular form using orthogonal
c                      transformations. usually followed by qzit,
c                      qzval, qz
c
c qzit     -      -    the second step of the qz algorithm for
c                      generalized eigenproblems.  accepts an upper
c                      hessenberg and an upper triangular matrix
c                      and reduces the former to quasi-triangular
c                      form while preserving the form of the latter.
c                      usually preceded by qzhes and followed by qzval
c                      and qzvec.
c
c figi     -      -    transforms certain real non-symmetric
c                      tridiagonal matrix to symmetric tridiagonal
c                      matrix.
c
c figi2    -      -    transforms certain real non-symmetric
c                      tridiagonal matrix to symmetric tridiagonal
c                      matrix.
c
c reduc    -      -    reduces generalized symmetric eigenproblem
c                      a*x=(lambda)*b*x, to standard symmetric
c                      eigenproblem using cholesky factorization.
c
c reduc2   -      -    reduces certain generalized symmetric
c                      eigenproblems standard symmetric eigenproblem,
c                      using cholesky factorization.
c
c   -      -    comlr  computes eigenvalues of a complex upper
c                      hessenberg matrix using the modified lr method.
c
c   -      -    comlr2 computes eigenvalues and eigenvectors of
c                      complex upper hessenberg matrix using
c                      modified lr method.
c
c hqr      -    comqr  computes eigenvalues of a real (complex)
c                      upper hessenberg matrix using the qr method.
c
c hqr2     -    comqr2 computes eigenvalues and eigenvectors of
c                      real (complex) upper hessenberg matrix
c                      using qr method.
c
c invit    -    cinvit computes eigenvectors of real (complex)
c                      hessenberg matrix associated with specified
c                      eigenvalues by inverse iteration.
c
c qzval    -      -    the third step of the qz algorithm for
c                      generalized eigenproblems.  accepts a pair
c                      of real matrices, one quasi-triangular form
c                      and the other in upper triangular form and
c                      computes the eigenvalues of the associated
c                      eigenproblem.  usually preceded by qzhes,
c                      qzit, and followed by qzvec.
c
c bandv    -      -    forms eigenvectors of real symmetric band
c                      matrix associated with a set of ordered
c                      approximate eigenvalue by inverse iteration.
c
c qzvec    -      -    the optional fourth step of the qz algorithm
c                      for generalized eigenproblems.  accepts
c                      a matrix in quasi-triangular form and another
c                      in upper triangular and computes the
c                      eigenvectors of the triangular problem
c                      and transforms them back to the original
c                      coordinates usually preceded by qzhes, qzit,
c                      qzval.
c
c tinvit   -      -    eigenvectors of symmetric tridiagonal
c                      matrix corresponding to some specified
c                      eigenvalues, using inverse iteration.
c
c bakvec   -      -    forms eigenvectors of certain real
c                      non-symmetric tridiagonal matrix from
c                      symmetric tridiagonal matrix output from figi.
c
c balbak   -    cbabk2 forms eigenvectors of real (complex) general
c                      matrix from eigenvectors of matrix output
c                      from balanc (cbal).
c
c elmbak   -    combak forms eigenvectors of real (complex) general
c                      matrix from eigenvectors of upper hessenberg
c                      matrix output from elmhes (comhes).
c
c eltran   -      -    accumulates the stabilized elementary
c                      similarity transformations used in the
c                      reduction of a real general matrix to upper
c                      hessenberg form by elmhes.
c
c htrib3   -      -    computes eigenvectors of complex hermitian
c                      matrix from eigenvectors of real symmetric
c                      tridiagonal matrix output from htrid3.
c
c htribk   -      -    forms eigenvectors of complex hermitian
c                      matrix from eigenvectors of real symmetric
c                      tridiagonal matrix output from htridi.
c
c ortbak   -    cortb  forms eigenvectors of general real (complex)
c                      matrix from eigenvectors of upper hessenberg
c                      matrix output from orthes (corth).
c
c ortran   -      -    accumulates orthogonal similarity
c                      transformations in reduction of real general
c                      matrix by orthes.
c
c rebak    -      -    forms eigenvectors of generalized symmetric
c                      eigensystem from eigenvectors of derived
c                      matrix output from reduc or reduc2.
c
c rebakb   -      -    forms eigenvectors of generalized symmetric
c                      eigensystem from eigenvectors of derived
c                      matrix output from reduc2
c
c trbak1   -      -    forms the eigenvectors of real symmetric
c                      matrix from eigenvectors of symmetric
c                      tridiagonal matrix formed by tred1.
c
c trbak3   -      -    forms eigenvectors of real symmetric matrix
c                      from the eigenvectors of symmetric tridiagonal
c                      matrix formed by tred3.
c
c minfit   -      -    compute singular value decomposition
c                      of rectangular matrix and solve related
c                      linear least squares problem.
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  (none)
c***revision history  (yymmdd)
c   811101  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900723  purpose section revised.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  eisdoc
c***first executable statement  eisdoc
      return
      end
