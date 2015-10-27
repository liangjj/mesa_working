*deck rt
      subroutine rt (nm, n, a, w, matz, z, fv1, ierr)
c***begin prologue  rt
c***purpose  compute the eigenvalues and eigenvectors of a special real
c            tridiagonal matrix.
c***library   slatec (eispack)
c***category  d4a5
c***type      single precision (rt-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine calls the recommended sequence of subroutines
c     from the eigensystem subroutine package (eispack) to find the
c     eigenvalues and eigenvectors (if desired) of a special real
c     tridiagonal matrix.  the property of the matrix required for use
c     of this subroutine is that the products of pairs of corresponding
c     off-diagonal elements be all non-negative.  if eigenvectors are
c     desired, no product can be zero unless both factors are zero.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameter, a and z, as declared in the calling
c          program dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix a.  n is an integer variable.
c          n must be less than or equal to nm.
c
c        a contains the special real tridiagonal matrix in its first
c          three columns.  the subdiagonal elements are stored in the
c          last n-1 positions of the first column, the diagonal elements
c          in the second column, and the superdiagonal elements in the
c          first n-1 positions of the third column.  elements a(1,1) and
c          a(n,3) are arbitrary.  a is a two-dimensional real array,
c          dimensioned a(nm,3).
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
c        z contains the eigenvectors if matz is not zero.  the eigen-
c          vectors are not normalized.  z is a two-dimensional real
c          array, dimensioned z(nm,n).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          10*n       if n is greater than nm,
c          n+j        if a(j,1)*a(j-1,3) is negative,
c          2*n+j      if the product is zero with one factor non-zero,
c                     and matz is non-zero;
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c                     the eigenvalues and eigenvectors in the w and z
c                     arrays should be correct for indices
c                     1, 2, ..., ierr-1.
c
c        fv1 is a one-dimensional real array used for temporary storage,
c          dimensioned fv1(n).
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  figi, figi2, imtql1, imtql2
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  rt
c
      integer n,nm,ierr,matz
      real a(nm,3),w(*),z(nm,*),fv1(*)
c
c***first executable statement  rt
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  figi(nm,n,a,w,fv1,fv1,ierr)
      if (ierr .gt. 0) go to 50
      call  imtql1(n,w,fv1,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  figi2(nm,n,a,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  imtql2(nm,n,w,fv1,z,ierr)
   50 return
      end
