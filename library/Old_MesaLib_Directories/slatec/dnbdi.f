*deck dnbdi
      subroutine dnbdi (abe, lda, n, ml, mu, ipvt, det)
c***begin prologue  dnbdi
c***purpose  compute the determinant of a band matrix using the factors
c            computed by dnbco or dnbfa.
c***library   slatec
c***category  d3a2
c***type      double precision (snbdi-s, dnbdi-d, cnbdi-c)
c***keywords  banded, determinant, linear equations, nonsymmetric
c***author  voorhees, e. a., (lanl)
c***description
c
c     dnbdi computes the determinant of a band matrix
c     using the factors computed by dnbco or dnbfa.
c     if the inverse is needed, use dnbsl  n  times.
c
c     on entry
c
c        abe     double precision(lda, nc)
c                the output from dnbco or dnbfa.
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
c                the pivot vector from dnbco or dnbfa.
c
c     on return
c
c        det     double precision(2)
c                determinant of original matrix.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. abs(det(1)) .lt. 10.0
c                or  det(1) = 0.0 .
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  (none)
c***revision history  (yymmdd)
c   800728  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dnbdi
      integer lda,n,ml,mu,ipvt(*)
      double precision abe(lda,*),det(2)
c
      double precision ten
      integer i
c***first executable statement  dnbdi
      det(1) = 1.0d0
      det(2) = 0.0d0
      ten = 10.0d0
      do 50 i = 1, n
         if (ipvt(i) .ne. i) det(1) = -det(1)
         det(1) = abe(i,ml+1)*det(1)
         if (det(1) .eq. 0.0d0) go to 60
   10    if (abs(det(1)) .ge. 1.0d0) go to 20
            det(1) = ten*det(1)
            det(2) = det(2) - 1.0d0
         go to 10
   20    continue
   30    if (abs(det(1)) .lt. ten) go to 40
            det(1) = det(1)/ten
            det(2) = det(2) + 1.0d0
         go to 30
   40    continue
   50 continue
   60 continue
      return
      end
