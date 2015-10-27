*deck cpofa
      subroutine cpofa (a, lda, n, info)
c***begin prologue  cpofa
c***purpose  factor a complex hermitian positive definite matrix.
c***library   slatec (linpack)
c***category  d2d1b
c***type      complex (spofa-s, dpofa-d, cpofa-c)
c***keywords  linear algebra, linpack, matrix factorization,
c             positive definite
c***author  moler, c. b., (u. of new mexico)
c***description
c
c     cpofa factors a complex hermitian positive definite matrix.
c
c     cpofa is usually called by cpoco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for cpoco) = (1 + 18/n)*(time for cpofa) .
c
c     on entry
c
c        a       complex(lda, n)
c                the hermitian matrix to be factored.  only the
c                diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix  r  so that  a =
c                ctrans(r)*r where  ctrans(r)  is the conjugate
c                transpose.  the strict lower triangle is unaltered.
c                if  info .ne. 0 , the factorization is not complete.
c
c        info    integer
c                = 0  for normal return.
c                = k  signals an error condition.  the leading minor
c                     of order  k  is not positive definite.
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  cdotc
c***revision history  (yymmdd)
c   780814  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cpofa
      integer lda,n,info
      complex a(lda,*)
c
      complex cdotc,t
      real s
      integer j,jm1,k
c***first executable statement  cpofa
         do 30 j = 1, n
            info = j
            s = 0.0e0
            jm1 = j - 1
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               t = a(k,j) - cdotc(k-1,a(1,k),1,a(1,j),1)
               t = t/a(k,k)
               a(k,j) = t
               s = s + real(t*conjg(t))
   10       continue
   20       continue
            s = real(a(j,j)) - s
            if (s .le. 0.0e0 .or. aimag(a(j,j)) .ne. 0.0e0) go to 40
            a(j,j) = cmplx(sqrt(s),0.0e0)
   30    continue
         info = 0
   40 continue
      return
      end
