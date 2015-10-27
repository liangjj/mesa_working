*deck cppsl
      subroutine cppsl (ap, n, b)
c***begin prologue  cppsl
c***purpose  solve the complex hermitian positive definite system using
c            the factors computed by cppco or cppfa.
c***library   slatec (linpack)
c***category  d2d1b
c***type      complex (sppsl-s, dppsl-d, cppsl-c)
c***keywords  linear algebra, linpack, matrix, packed,
c             positive definite, solve
c***author  moler, c. b., (u. of new mexico)
c***description
c
c     cppsl solves the complex hermitian positive definite system
c     a * x = b
c     using the factors computed by cppco or cppfa.
c
c     on entry
c
c        ap      complex (n*(n+1)/2)
c                the output from cppco or cppfa.
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
c           call cppco(ap,n,rcond,z,info)
c           if (rcond is too small .or. info .ne. 0) go to ...
c           do 10 j = 1, p
c              call cppsl(ap,n,c(1,j))
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
c***end prologue  cppsl
      integer n
      complex ap(*),b(*)
c
      complex cdotc,t
      integer k,kb,kk
c***first executable statement  cppsl
      kk = 0
      do 10 k = 1, n
         t = cdotc(k-1,ap(kk+1),1,b(1),1)
         kk = kk + k
         b(k) = (b(k) - t)/ap(kk)
   10 continue
      do 20 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/ap(kk)
         kk = kk - k
         t = -b(k)
         call caxpy(k-1,t,ap(kk+1),1,b(1),1)
   20 continue
      return
      end
