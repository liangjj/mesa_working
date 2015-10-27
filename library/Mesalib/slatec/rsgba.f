*deck rsgba
      subroutine rsgba (nm, n, a, b, w, matz, z, fv1, fv2, ierr)
c***begin prologue  rsgba
c***purpose  compute the eigenvalues and, optionally, the eigenvectors
c            of a symmetric generalized eigenproblem.
c***library   slatec (eispack)
c***category  d4b1
c***type      single precision (rsgba-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     for the real symmetric generalized eigenproblem  bax = (lambda)x.
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
c        a contains a real symmetric matrix.  a is a two-dimensional
c          real array, dimensioned a(nm,n).
c
c        b contains a positive definite real symmetric matrix.  b is a
c          two-dimensional real array, dimensioned b(nm,n).
c
c        matz is an integer variable set equal to zero if only
c          eigenvalues are desired.  otherwise, it is set to any
c          non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w contains the eigenvalues in ascending order.  w is a
c          one-dimensional real array, dimensioned w(n).
c
c        z contains the eigenvectors if matz is not zero.  z is a
c          two-dimensional real array, dimensioned z(nm,n).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          10*n       if n is greater than nm,
c          7*n+1      if b is not positive definite,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c                     the eigenvalues should be correct for indices
c                     1, 2, ..., ierr-1, but no eigenvectors are
c                     computed.
c
c        fv1 and fv2 are one-dimensional real arrays used for temporary
c          storage, dimensioned fv1(n) and fv2(n).
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  rebakb, reduc2, tql2, tqlrat, tred1, tred2
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  rsgba
c
      integer n,nm,ierr,matz
      real a(nm,*),b(nm,*),w(*),z(nm,*),fv1(*),fv2(*)
c
c***first executable statement  rsgba
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  reduc2(nm,n,a,b,fv2,ierr)
      if (ierr .ne. 0) go to 50
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  rebakb(nm,n,b,fv2,n,z)
   50 return
      end
