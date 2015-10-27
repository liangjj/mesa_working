*deck rsb
      subroutine rsb (nm, n, mb, a, w, matz, z, fv1, fv2, ierr)
c***begin prologue  rsb
c***purpose  compute the eigenvalues and, optionally, the eigenvectors
c            of a symmetric band matrix.
c***library   slatec (eispack)
c***category  d4a6
c***type      single precision (rsb-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric band matrix.
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
c        mb is the half band width of the matrix, defined as the
c          number of adjacent diagonals, including the principal
c          diagonal, required to specify the non-zero portion of the
c          lower triangle of the matrix.  mb must be less than or
c          equal to n.  mb is an integer variable.
c
c        a contains the lower triangle of the real symmetric band
c          matrix.  its lowest subdiagonal is stored in the last
c          n+1-mb  positions of the first column, its next subdiagonal
c          in the last  n+2-mb  positions of the second column, further
c          subdiagonals similarly, and finally its principal diagonal
c          in the  n  positions of the last column.  contents of storage
c          locations not part of the matrix are arbitrary.  a is a
c          two-dimensional real array, dimensioned a(nm,mb).
c
c        matz is an integer variable set equal to zero if only
c          eigenvalues are desired.  otherwise, it is set to any
c          non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        a has been destroyed.
c
c        w contains the eigenvalues in ascending order.  w is a one-
c          dimensional real array, dimensioned w(n).
c
c        z contains the eigenvectors if matz is not zero.  the
c          eigenvectors are orthonormal.  z is a two-dimensional
c          real array, dimensioned z(nm,n).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          10*n       if n is greater than nm,
c          12*n       if mb is either non-positive or greater than n,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c                     the eigenvalues and eigenvectors, if requested,
c                     should be correct for indices 1, 2, ..., ierr-1.
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
c***routines called  bandr, tql2, tqlrat
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  rsb
c
      integer n,mb,nm,ierr,matz
      real a(nm,*),w(*),z(nm,*),fv1(*),fv2(*)
      logical tf
c
c***first executable statement  rsb
      if (n .le. nm) go to 5
      ierr = 10 * n
      go to 50
    5 if (mb .gt. 0) go to 10
      ierr = 12 * n
      go to 50
   10 if (mb .le. n) go to 15
      ierr = 12 * n
      go to 50
c
   15 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      tf = .false.
      call  bandr(nm,n,mb,a,w,fv1,fv2,tf,z)
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 tf = .true.
      call  bandr(nm,n,mb,a,w,fv1,fv1,tf,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end
