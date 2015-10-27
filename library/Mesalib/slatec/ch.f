*deck ch
      subroutine ch (nm, n, ar, ai, w, matz, zr, zi, fv1, fv2, fm1,
     +   ierr)
c***begin prologue  ch
c***purpose  compute the eigenvalues and, optionally, the eigenvectors
c            of a complex hermitian matrix.
c***library   slatec (eispack)
c***category  d4a3
c***type      complex (rs-s, ch-c)
c***keywords  eigenvalues, eigenvectors, eispack
c***author  smith, b. t., et al.
c***description
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a complex hermitian matrix.
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
c          of the complex hermitian matrix.  ar and ai are
c          two-dimensional real arrays, dimensioned ar(nm,n)
c          and ai(nm,n).
c
c        matz is an integer variable set equal to zero if only
c          eigenvalues are desired.  otherwise, it is set to any
c          non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w contains the eigenvalues in ascending order.
c          w is a one-dimensional real array, dimensioned w(n).
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
c                     1, 2, ..., ierr-1, but no eigenvectors are
c                     computed.
c
c        fv1 and fv2 are one-dimensional real arrays used for
c          temporary storage, dimensioned fv1(n) and fv2(n).
c
c        fm1 is a two-dimensional real array used for temporary
c          storage, dimensioned fm1(2,n).
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c     ------------------------------------------------------------------
c
c***references  b. t. smith, j. m. boyle, j. j. dongarra, b. s. garbow,
c                 y. ikebe, v. c. klema and c. b. moler, matrix eigen-
c                 system routines - eispack guide, springer-verlag,
c                 1976.
c***routines called  htribk, htridi, tql2, tqlrat
c***revision history  (yymmdd)
c   760101  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  ch
c
      integer i,j,n,nm,ierr,matz
      real ar(nm,*),ai(nm,*),w(*),zr(nm,*),zi(nm,*)
      real fv1(*),fv2(*),fm1(2,*)
c
c***first executable statement  ch
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  htridi(nm,n,ar,ai,w,fv1,fv2,fm1)
      if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tqlrat(n,w,fv2,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 do 40 i = 1, n
c
         do 30 j = 1, n
            zr(j,i) = 0.0e0
   30    continue
c
         zr(i,i) = 1.0e0
   40 continue
c
      call  tql2(nm,n,w,fv1,zr,ierr)
      if (ierr .ne. 0) go to 50
      call  htribk(nm,n,ar,ai,fm1,n,zr,zi)
   50 return
      end
