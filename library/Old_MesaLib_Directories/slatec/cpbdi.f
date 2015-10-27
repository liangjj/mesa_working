*deck cpbdi
      subroutine cpbdi (abd, lda, n, m, det)
c***begin prologue  cpbdi
c***purpose  compute the determinant of a complex hermitian positive
c            definite band matrix using the factors computed by cpbco or
c            cpbfa.
c***library   slatec (linpack)
c***category  d3d2
c***type      complex (spbdi-s, dpbdi-d, cpbdi-c)
c***keywords  banded, determinant, inverse, linear algebra, linpack,
c             matrix, positive definite
c***author  moler, c. b., (u. of new mexico)
c***description
c
c     cpbdi computes the determinant
c     of a complex hermitian positive definite band matrix
c     using the factors computed by cpbco or cpbfa.
c     if the inverse is needed, use cpbsl  n  times.
c
c     on entry
c
c        abd     complex(lda, n)
c                the output from cpbco or cpbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the matrix  a .
c
c        m       integer
c                the number of diagonals above the main diagonal.
c
c     on return
c
c        det     real(2)
c                determinant of original matrix in the form
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. det(1) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  (none)
c***revision history  (yymmdd)
c   780814  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cpbdi
      integer lda,n,m
      complex abd(lda,*)
      real det(2)
c
      real s
      integer i
c***first executable statement  cpbdi
c
c     compute determinant
c
      det(1) = 1.0e0
      det(2) = 0.0e0
      s = 10.0e0
      do 50 i = 1, n
         det(1) = real(abd(m+1,i))**2*det(1)
         if (det(1) .eq. 0.0e0) go to 60
   10    if (det(1) .ge. 1.0e0) go to 20
            det(1) = s*det(1)
            det(2) = det(2) - 1.0e0
         go to 10
   20    continue
   30    if (det(1) .lt. s) go to 40
            det(1) = det(1)/s
            det(2) = det(2) + 1.0e0
         go to 30
   40    continue
   50 continue
   60 continue
      return
      end
