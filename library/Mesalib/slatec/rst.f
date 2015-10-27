*deck rst
      subroutine rst (nm, n, w, e, matz, z, ierr)
c***begin prologue  rst
c***purpose  compute the eigenvalues and, optionally, the eigenvectors
c            of a real symmetric tridiagonal matrix.
c***library   slatec (eispack)
c***category  d4a5
c***type      single precision (rst-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric tridiagonal matrix.
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
c        w contains the diagonal elements of the real symmetric
c          tridiagonal matrix.  w is a one-dimensional real array,
c          dimensioned w(n).
c
c        e contains the subdiagonal elements of the matrix in its last
c          n-1 positions.  e(1) is arbitrary.  e is a one-dimensional
c          real array, dimensioned e(n).
c
c        matz is an integer variable set equal to zero if only
c          eigenvalues are desired.  otherwise, it is set to any
c          non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w contains the eigenvalues in ascending order.
c
c        z contains the eigenvectors if matz is not zero.  the eigen-
c          vectors are orthonormal.  z is a two-dimensional real array,
c          dimensioned z(nm,n).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          10*n       if n is greater than nm,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c                     the eigenvalues and eigenvectors in the w and z
c                     arrays should be correct for indices
c                     1, 2, ..., ierr-1.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  imtql1, imtql2
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  rst
c
      integer i,j,n,nm,ierr,matz
      real w(*),e(*),z(nm,*)
c
c***first executable statement  rst
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  imtql1(n,w,e,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
c
         do 30 j = 1, n
            z(j,i) = 0.0e0
   30    continue
c
         z(i,i) = 1.0e0
   40 continue
c
      call  imtql2(nm,n,w,e,z,ierr)
   50 return
      end
