*deck rsp
      subroutine rsp (nm, n, nv, a, w, matz, z, fv1, fv2, ierr)
c***begin prologue  rsp
c***purpose  compute the eigenvalues and, optionally, the eigenvectors
c            of a real symmetric matrix packed into a one dimensional
c            array.
c***library   slatec (eispack)
c***category  d4a1
c***type      single precision (rsp-s)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric packed matrix.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameter, z, as declared in the calling program
c          dimension statement.  nm is an integer variable.
c
c        n is the order of the matrix a.  n is an integer variable.
c          n must be less than or equal to nm.
c
c        nv is an integer variable set equal to the dimension of the
c          array a as specified in the calling program.  nv must not
c          be less than  n*(n+1)/2.
c
c        a contains the lower triangle, stored row-wise, of the real
c          symmetric packed matrix.  a is a one-dimensional real
c          array, dimensioned a(nv).
c
c        matz is an integer variable set equal to zero if only
c          eigenvalues are desired.  otherwise, it is set to any
c          non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        a has been destroyed.
c
c        w contains the eigenvalues in ascending order.  w is a
c          one-dimensional real array, dimensioned w(n).
c
c        z contains the eigenvectors if matz is not zero.  the eigen-
c          vectors are orthonormal.  z is a two-dimensional real array,
c          dimensioned z(nm,n).
c
c        ierr is an integer flag set to
c          zero       for normal return,
c          10*n       if n is greater than nm,
c          20*n       if nv is less than n*(n+1)/2,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c                     the eigenvalues and eigenvectors in the w and z
c                     arrays should be correct for indices
c                     1, 2, ..., ierr-1.
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
c***routines called  tql2, tqlrat, trbak3, tred3
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  rsp
c
      integer i,j,n,nm,nv,ierr,matz
      real a(*),w(*),z(nm,*),fv1(*),fv2(*)
c
c***first executable statement  rsp
      if (n .le. nm) go to 5
      ierr = 10 * n
      go to 50
    5 if (nv .ge. (n * (n + 1)) / 2) go to 10
      ierr = 20 * n
      go to 50
c
   10 call  tred3(n,nv,a,w,fv1,fv2)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tqlrat(n,w,fv2,ierr)
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
      call  tql2(nm,n,w,fv1,z,ierr)
      if (ierr .ne. 0) go to 50
      call  trbak3(nm,n,nv,a,n,z)
   50 return
      end
