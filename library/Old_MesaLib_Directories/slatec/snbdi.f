*deck snbdi
      subroutine snbdi (abe, lda, n, ml, mu, ipvt, det)
c***begin prologue  snbdi
c***purpose  compute the determinant of a band matrix using the factors
c            computed by snbco or snbfa.
c***library   slatec
c***category  d3a2
c***type      single precision (snbdi-s, dnbdi-d, cnbdi-c)
c***keywords  banded, determinant, linear equations, nonsymmetric
c***author  voorhees, e. a., (lanl)
c***description
c
c     snbdi computes the determinant of a band matrix
c     using the factors computed by snbco or snbfa.
c     if the inverse is needed, use snbsl  n  times.
c
c     on entry
c
c        abe     real(lda, nc)
c                the output from snbco or snbfa.
c                nc must be .ge. 2*ml+mu+1 .
c
c        lda     integer
c                the leading dimension of the array  abe .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from snbco or snbfa.
c
c     on return
c
c        det     real(2)
c                determinant of original matrix.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. abs(det(1)) .lt. 10.0
c                or  det(1) = 0.0 .
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  (none)
c***revision history  (yymmdd)
c   800725  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  snbdi
      integer lda,n,ml,mu,ipvt(*)
      real abe(lda,*),det(2)
c
      real ten
      integer i
c***first executable statement  snbdi
      det(1) = 1.0e0
      det(2) = 0.0e0
      ten = 10.0e0
      do 50 i = 1, n
         if (ipvt(i) .ne. i) det(1) = -det(1)
         det(1) = abe(i,ml+1)*det(1)
         if (det(1) .eq. 0.0e0) go to 60
   10    if (abs(det(1)) .ge. 1.0e0) go to 20
            det(1) = ten*det(1)
            det(2) = det(2) - 1.0e0
         go to 10
   20    continue
   30    if (abs(det(1)) .lt. ten) go to 40
            det(1) = det(1)/ten
            det(2) = det(2) + 1.0e0
         go to 30
   40    continue
   50 continue
   60 continue
      return
      end
