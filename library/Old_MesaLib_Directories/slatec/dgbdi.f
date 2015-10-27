*deck dgbdi
      subroutine dgbdi (abd, lda, n, ml, mu, ipvt, det)
c***begin prologue  dgbdi
c***purpose  compute the determinant of a band matrix using the factors
c            computed by dgbco or dgbfa.
c***library   slatec (linpack)
c***category  d3a2
c***type      double precision (sgbdi-s, dgbdi-d, cgbdi-c)
c***keywords  banded, determinant, inverse, linear algebra, linpack,
c             matrix
c***author  moler, c. b., (u. of new mexico)
c***description
c
c     dgbdi computes the determinant of a band matrix
c     using the factors computed by dgbco or dgbfa.
c     if the inverse is needed, use dgbsl  n  times.
c
c     on entry
c
c        abd     double precision(lda, n)
c                the output from dgbco or dgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
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
c                the pivot vector from dgbco or dgbfa.
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
c   780814  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dgbdi
      integer lda,n,ml,mu,ipvt(*)
      double precision abd(lda,*),det(2)
c
      double precision ten
      integer i,m
c***first executable statement  dgbdi
      m = ml + mu + 1
      det(1) = 1.0d0
      det(2) = 0.0d0
      ten = 10.0d0
      do 50 i = 1, n
         if (ipvt(i) .ne. i) det(1) = -det(1)
         det(1) = abd(m,i)*det(1)
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
