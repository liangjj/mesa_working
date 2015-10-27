*deck dpbsl
      subroutine dpbsl (abd, lda, n, m, b)
c***begin prologue  dpbsl
c***purpose  solve a real symmetric positive definite band system
c            using the factors computed by dpbco or dpbfa.
c***library   slatec (linpack)
c***category  d2b2
c***type      double precision (spbsl-s, dpbsl-d, cpbsl-c)
c***keywords  banded, linear algebra, linpack, matrix,
c             positive definite, solve
c***author  moler, c. b., (u. of new mexico)
c***description
c
c     dpbsl solves the double precision symmetric positive definite
c     band system  a*x = b
c     using the factors computed by dpbco or dpbfa.
c
c     on entry
c
c        abd     double precision(lda, n)
c                the output from dpbco or dpbfa.
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
c        b       double precision(n)
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
c        singularity, but it is usually caused by improper subroutine
c        arguments.  it will not occur if the subroutines are called
c        correctly, and  info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dpbco(abd,lda,n,rcond,z,info)
c           if (rcond is too small .or. info .ne. 0) go to ...
c           do 10 j = 1, p
c              call dpbsl(abd,lda,n,c(1,j))
c        10 continue
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  daxpy, ddot
c***revision history  (yymmdd)
c   780814  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dpbsl
      integer lda,n,m
      double precision abd(lda,*),b(*)
c
      double precision ddot,t
      integer k,kb,la,lb,lm
c
c     solve trans(r)*y = b
c
c***first executable statement  dpbsl
      do 10 k = 1, n
         lm = min(k-1,m)
         la = m + 1 - lm
         lb = k - lm
         t = ddot(lm,abd(la,k),1,b(lb),1)
         b(k) = (b(k) - t)/abd(m+1,k)
   10 continue
c
c     solve r*x = y
c
      do 20 kb = 1, n
         k = n + 1 - kb
         lm = min(k-1,m)
         la = m + 1 - lm
         lb = k - lm
         b(k) = b(k)/abd(m+1,k)
         t = -b(k)
         call daxpy(lm,t,abd(la,k),1,b(lb),1)
   20 continue
      return
      end