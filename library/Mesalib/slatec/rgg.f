*deck rgg
      subroutine rgg (nm, n, a, b, alfr, alfi, beta, matz, z, ierr)
c***begin prologue  rgg
c***purpose  compute the eigenvalues and eigenvectors for a real
c            generalized eigenproblem.
c***library   slatec (eispack)
c***category  d4b2
c***type      single precision (rgg-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     for the real general generalized eigenproblem  ax = (lambda)bx.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, a, b, and z, as declared in the calling
c          program dimension statement.  nm is an integer variable.
c
c        n is the order of the matrices a and b.  n is an integer
c          variable.  n must be less than or equal to nm.
c
c        a contains a real general matrix.  a is a two-dimensional
c          real array, dimensioned a(nm,n).
c
c        b contains a real general matrix.  b is a two-dimensional
c          real array, dimensioned b(nm,n).
c
c        matz is an integer variable set equal to zero if only
c          eigenvalues are desired.  otherwise, it is set to any
c          non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        a and b have been destroyed.
c
c        alfr and alfi contain the real and imaginary parts,
c          respectively, of the numerators of the eigenvalues.
c          alfr and alfi are one-dimensional real arrays,
c          dimensioned alfr(n) and alfi(n).
c
c        beta contains the denominators of the eigenvalues,
c          which are thus given by the ratios  (alfr+i*alfi)/beta.
c          complex conjugate pairs of eigenvalues appear consecutively
c          with the eigenvalue having the positive imaginary part first.
c          beta is a one-dimensional real array, dimensioned beta(n).
c
c        z contains the real and imaginary parts of the eigenvectors
c          if matz is not zero.  if the j-th eigenvalue is real, the
c          j-th column of  z  contains its eigenvector.  if the j-th
c          eigenvalue is complex with positive imaginary part, the
c          j-th and (j+1)-th columns of  z  contain the real and
c          imaginary parts of its eigenvector.  the conjugate of this
c          vector is the eigenvector for the conjugate eigenvalue.
c          z is a two-dimensional real array, dimensioned z(nm,n).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          10*n       if n is greater than nm,
c          j          if the j-th eigenvalue has not been
c                     determined after a total of 30*n iterations.
c                     the eigenvalues should be correct for indices
c                     ierr+1, ierr+2, ..., n, but no eigenvectors are
c                     computed.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  qzhes, qzit, qzval, qzvec
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  rgg
c
      integer n,nm,ierr,matz
      real a(nm,*),b(nm,*),alfr(*),alfi(*),beta(*),z(nm,*)
      logical tf
c
c***first executable statement  rgg
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      tf = .false.
      call  qzhes(nm,n,a,b,tf,z)
      call  qzit(nm,n,a,b,0.0e0,tf,z,ierr)
      call  qzval(nm,n,a,b,alfr,alfi,beta,tf,z)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 tf = .true.
      call  qzhes(nm,n,a,b,tf,z)
      call  qzit(nm,n,a,b,0.0e0,tf,z,ierr)
      call  qzval(nm,n,a,b,alfr,alfi,beta,tf,z)
      if (ierr .ne. 0) go to 50
      call  qzvec(nm,n,a,b,alfr,alfi,beta,z)
   50 return
      end
