*deck cg
      subroutine cg (nm, n, ar, ai, wr, wi, matz, zr, zi, fv1, fv2, fv3,
     +   ierr)
c***begin prologue  cg
c***purpose  compute the eigenvalues and, optionally, the eigenvectors
c            of a complex general matrix.
c***library   slatec (eispack)
c***category  d4a4
c***type      complex (rg-s, cg-c)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex general matrix.
c
c     on input
c
c        nm must be set to the row dimension of the two-dimensional
c          array parameters, ar, ai, zr and zi, as declared in the
c          calling program dimension statement.  nm is an integer
c          variable.
c
c        n is the order of the matrix a=(ar,ai).  n is an integer
c          variable.  n must be less than or equal to nm.
c
c        ar and ai contain the real and imaginary parts, respectively,
c          of the complex general matrix.  ar and ai are two-dimensional
c          real arrays, dimensioned ar(nm,n) and ai(nm,n).
c
c        matz is an integer variable set equal to zero if only
c          eigenvalues are desired.  otherwise, it is set to any
c          non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        wr and wi contain the real and imaginary parts, respectively,
c          of the eigenvalues.  wr and wi are one-dimensional real
c          arrays, dimensioned wr(n) and wi(n).
c
c        zr and zi contain the real and imaginary parts, respectively,
c          of the eigenvectors if matz is not zero.  zr and zi are
c          two-dimensional real arrays, dimensioned zr(nm,n) and
c          zi(nm,n).
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
c        fv1, fv2, and fv3 are one-dimensional real arrays used for
c          temporary storage, dimensioned fv1(n), fv2(n), and fv3(n).
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  cbabk2, cbal, comqr, comqr2, corth
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cg
c
      integer n,nm,is1,is2,ierr,matz
      real ar(nm,*),ai(nm,*),wr(*),wi(*),zr(nm,*),zi(nm,*)
      real fv1(*),fv2(*),fv3(*)
c
c***first executable statement  cg
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  cbal(nm,n,ar,ai,is1,is2,fv1)
      call  corth(nm,n,is1,is2,ar,ai,fv2,fv3)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  comqr(nm,n,is1,is2,ar,ai,wr,wi,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  comqr2(nm,n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi,ierr)
      if (ierr .ne. 0) go to 50
      call  cbabk2(nm,n,is1,is2,fv1,n,zr,zi)
   50 return
      end
