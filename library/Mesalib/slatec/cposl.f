*deck cposl
      subroutine cposl (a, lda, n, b)
c***begin prologue  cposl
c***purpose  solve the complex hermitian positive definite linear system
c            using the factors computed by cpoco or cpofa.
c***library   slatec (linpack)
c***category  d2d1b
c***type      complex (sposl-s, dposl-d, cposl-c)
c***keywords  linear algebra, linpack, matrix, positive definite, solve
c***author  moler, c. b., (u. of new mexico)
c***description
c
c     cposl solves the complex hermitian positive definite system
c     a * x = b
c     using the factors computed by cpoco or cpofa.
c
c     on entry
c
c        a       complex(lda, n)
c                the output from cpoco or cpofa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        b       complex(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal.  technically this indicates
c        singularity but it is usually caused by improper subroutine
c        arguments.  it will not occur if the subroutines are called
c        correctly and  info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call cpoco(a,lda,n,rcond,z,info)
c           if (rcond is too small .or. info .ne. 0) go to ...
c           do 10 j = 1, p
c              call cposl(a,lda,n,c(1,j))
c        10 continue
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  caxpy, cdotc
c***revision history  (yymmdd)
c   780814  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cposl
      integer lda,n
      complex a(lda,*),b(*)
c
      complex cdotc,t
      integer k,kb
c
c     solve ctrans(r)*y = b
c
c***first executable statement  cposl
      do 10 k = 1, n
         t = cdotc(k-1,a(1,k),1,b(1),1)
         b(k) = (b(k) - t)/a(k,k)
   10 continue
c
c     solve r*x = y
c
      do 20 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/a(k,k)
         t = -b(k)
         call caxpy(k-1,t,a(1,k),1,b(1),1)
   20 continue
      return
      end
