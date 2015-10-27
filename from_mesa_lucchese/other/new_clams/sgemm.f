*deck %W% %G%
      subroutine sgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb,
     $                   beta, c, ldc )
c     .. scalar arguments ..
      character*1        transa, transb
      integer            m, n, k, lda, ldb, ldc
      real*8             alpha, beta
c     .. array arguments ..
      real*8             a( lda, * ), b( ldb, * ), c( ldc, * )
c     ..
c
c  purpose
c  =======
c
c  sgemm  performs one of the matrix-matrix operations
c
c     c := alpha*op( a )*op( b ) + beta*c,
c
c  where  op( x ) is one of
c
c     op( x ) = x   or   op( x ) = x',
c
c  alpha and beta are scalars, and a, b and c are matrices, with op( a )
c  an m by k matrix,  op( b )  a  k by n matrix and  c an m by n matrix.
c
c  parameters
c  ==========
c
c  transa - character*1.
c           on entry, transa specifies the form of op( a ) to be used in
c           the matrix multiplication as follows:
c
c              transa = 'n' or 'n',  op( a ) = a.
c
c              transa = 't' or 't',  op( a ) = a'.
c
c              transa = 'c' or 'c',  op( a ) = a'.
c
c           unchanged on exit.
c
c  transb - character*1.
c           on entry, transb specifies the form of op( b ) to be used in
c           the matrix multiplication as follows:
c
c              transb = 'n' or 'n',  op( b ) = b.
c
c              transb = 't' or 't',  op( b ) = b'.
c
c              transb = 'c' or 'c',  op( b ) = b'.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry,  m  specifies  the number  of rows  of the  matrix
c           op( a )  and of the  matrix  c.  m  must  be at least  zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry,  n  specifies the number  of columns of the matrix
c           op( b ) and the number of columns of the matrix c. n must be
c           at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry,  k  specifies  the number of columns of the matrix
c           op( a ) and the number of rows of the matrix op( b ). k must
c           be at least  zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, ka ), where ka is
c           k  when  transa = 'n' or 'n',  and is  m  otherwise.
c           before entry with  transa = 'n' or 'n',  the leading  m by k
c           part of the array  a  must contain the matrix  a,  otherwise
c           the leading  k by m  part of the array  a  must contain  the
c           matrix a.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. when  transa = 'n' or 'n' then
c           lda must be at least  max( 1, m ), otherwise  lda must be at
c           least  max( 1, k ).
c           unchanged on exit.
c
c  b      - real             array of dimension ( ldb, kb ), where kb is
c           n  when  transb = 'n' or 'n',  and is  k  otherwise.
c           before entry with  transb = 'n' or 'n',  the leading  k by n
c           part of the array  b  must contain the matrix  b,  otherwise
c           the leading  n by k  part of the array  b  must contain  the
c           matrix b.
c           unchanged on exit.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in the calling (sub) program. when  transb = 'n' or 'n' then
c           ldb must be at least  max( 1, k ), otherwise  ldb must be at
c           least  max( 1, n ).
c           unchanged on exit.
c
c  beta   - real            .
c           on entry,  beta  specifies the scalar  beta.  when  beta  is
c           supplied as zero then c need not be set on input.
c           unchanged on exit.
c
c  c      - real             array of dimension ( ldc, n ).
c           before entry, the leading  m by n  part of the array  c must
c           contain the matrix  c,  except when  beta  is zero, in which
c           case c need not be set on entry.
c           on exit, the array  c  is overwritten by the  m by n  matrix
c           ( alpha*op( a )*op( b ) + beta*c ).
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
      call dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb,
     $                   beta, c, ldc )
      return
c
c     end of sgemm .
c
      end
