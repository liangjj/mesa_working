*deck rg
      subroutine rg (nm, n, a, wr, wi, matz, z, iv1, fv1, ierr)
c***begin prologue  rg
c***purpose  compute the eigenvalues and, optionally, the eigenvectors
c            of a real general matrix.
c***library   slatec (eispack)
c***category  d4a2
c***type      single precision (rg-s, cg-c)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real general matrix.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, a and z, as declared in the calling
c          program dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix a.  n is an integer variable.
c          n must be less than or equal to nm.
c
c        a contains the real general matrix.  a is a two-dimensional
c          real array, dimensioned a(nm,n).
c
c        matz is an integer variable set equal to zero if only
c          eigenvalues are desired.  otherwise, it is set to any
c          non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        a has been destroyed.
c
c        wr and wi contain the real and imaginary parts, respectively,
c          of the eigenvalues.  the eigenvalues are unordered except
c          that complex conjugate pairs of eigenvalues appear consecu-
c          tively with the eigenvalue having the positive imaginary part
c          first.  if an error exit is made, the eigenvalues should be
c          correct for indices ierr+1, ierr+2, ..., n.  wr and wi are
c          one-dimensional real arrays, dimensioned wr(n) and wi(n).
c
c        z contains the real and imaginary parts of the eigenvectors
c          if matz is not zero.  if the j-th eigenvalue is real, the
c          j-th column of z contains its eigenvector.  if the j-th
c          eigenvalue is complex with positive imaginary part, the
c          j-th and (j+1)-th columns of z contain the real and
c          imaginary parts of its eigenvector.  the conjugate of this
c          vector is the eigenvector for the conjugate eigenvalue.
c          z is a two-dimensional real array, dimensioned z(nm,n).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          10*n       if n is greater than nm,
c          j          if the j-th eigenvalue has not been
c                     determined after a total of 30 iterations.
c                     the eigenvalues should be correct for indices
c                     ierr+1, ierr+2, ..., n, but no eigenvectors are
c                     computed.
c
c        iv1 and fv1 are one-dimensional temporary storage arrays of
c          dimension n.  iv1 is of type integer and fv1 of type real.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  balanc, balbak, elmhes, eltran, hqr, hqr2
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c   921103  corrected description of iv1.  (dwl, fnf and wrb)
c***end prologue  rg
c
      integer n,nm,is1,is2,ierr,matz
      real a(nm,*),wr(*),wi(*),z(nm,*),fv1(*)
      integer iv1(*)
c
c***first executable statement  rg
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  balanc(nm,n,a,is1,is2,fv1)
      call  elmhes(nm,n,is1,is2,a,iv1)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  hqr(nm,n,is1,is2,a,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  eltran(nm,n,is1,is2,a,iv1,z)
      call  hqr2(nm,n,is1,is2,a,wr,wi,z,ierr)
      if (ierr .ne. 0) go to 50
      call  balbak(nm,n,is1,is2,fv1,n,z)
   50 return
      end
