*deck cnbdi
      subroutine cnbdi (abe, lda, n, ml, mu, ipvt, det)
c***begin prologue  cnbdi
c***purpose  compute the determinant of a band matrix using the factors
c            computed by cnbco or cnbfa.
c***library   slatec
c***category  d3c2
c***type      complex (snbdi-s, dnbdi-d, cnbdi-c)
c***keywords  banded, determinant, linear equations, nonsymmetric
c***author  voorhees, e. a., (lanl)
c***description
c
c     cnbdi computes the determinant of a band matrix
c     using the factors computed by cnbco or cnbfa.
c     if the inverse is needed, use cnbsl  n  times.
c
c     on entry
c
c        abe     complex(lda, nc)
c                the output from cnbco or cnbfa.
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
c                the pivot vector from cnbco or cnbfa.
c
c     on return
c
c        det     complex(2)
c                determinant of original matrix.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. cabs1(det(1)) .lt. 10.0
c                or  det(1) = 0.0 .
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  (none)
c***revision history  (yymmdd)
c   800730  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cnbdi
      integer lda,n,ml,mu,ipvt(*)
      complex abe(lda,*),det(2)
c
      real ten
      integer i
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c***first executable statement  cnbdi
      det(1) = (1.0e0,0.0e0)
      det(2) = (0.0e0,0.0e0)
      ten = 10.0e0
      do 50 i = 1, n
         if (ipvt(i) .ne. i) det(1) = -det(1)
         det(1) = abe(i,ml+1)*det(1)
         if (cabs1(det(1)) .eq. 0.0e0) go to 60
   10    if (cabs1(det(1)) .ge. 1.0e0) go to 20
            det(1) = cmplx(ten,0.0e0)*det(1)
            det(2) = det(2) - (1.0e0,0.0e0)
         go to 10
   20    continue
   30    if (cabs1(det(1)) .lt. ten) go to 40
            det(1) = det(1)/cmplx(ten,0.0e0)
            det(2) = det(2) + (1.0e0,0.0e0)
         go to 30
   40    continue
   50 continue
   60 continue
      return
      end
