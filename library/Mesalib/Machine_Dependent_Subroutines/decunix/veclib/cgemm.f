*deck @(#)cgemm.f	1.4 8/6/91
      subroutine cgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb,
     $                   beta, c, ldc )
*     .. scalar arguments ..
      character*1        transa, transb
      integer            m, n, k, lda, ldb, ldc
      complex *16           alpha, beta
*     .. array arguments ..
      complex *16        a( lda, * ), b( ldb, * ), c( ldc, * )
*     ..
*
*  purpose
*  =======
*
*  cgemm  performs one of the matrix-matrix operations
*
*     c := alpha*op( a )*op( b ) + beta*c,
*
*  where  op( x ) is one of
*
*     op( x ) = x   or   op( x ) = x'   or   op( x ) = conjg( x' ),
*
*  alpha and beta are scalars, and a, b and c are matrices, with op( a )
*  an m by k matrix,  op( b )  a  k by n matrix and  c an m by n matrix.
*
*  parameters
*  ==========
*
*  transa - character*1.
*           on entry, transa specifies the form of op( a ) to be used in
*           the matrix multiplication as follows:
*
*              transa = 'n' or 'n',  op( a ) = a.
*
*              transa = 't' or 't',  op( a ) = a'.
*
*              transa = 'c' or 'c',  op( a ) = conjg( a' ).
*
*           unchanged on exit.
*
*  transb - character*1.
*           on entry, transb specifies the form of op( b ) to be used in
*           the matrix multiplication as follows:
*
*              transb = 'n' or 'n',  op( b ) = b.
*
*              transb = 't' or 't',  op( b ) = b'.
*
*              transb = 'c' or 'c',  op( b ) = conjg( b' ).
*
*           unchanged on exit.
*
*  m      - integer.
*           on entry,  m  specifies  the number  of rows  of the  matrix
*           op( a )  and of the  matrix  c.  m  must  be at least  zero.
*           unchanged on exit.
*
*  n      - integer.
*           on entry,  n  specifies the number  of columns of the matrix
*           op( b ) and the number of columns of the matrix c. n must be
*           at least zero.
*           unchanged on exit.
*
*  k      - integer.
*           on entry,  k  specifies  the number of columns of the matrix
*           op( a ) and the number of rows of the matrix op( b ). k must
*           be at least  zero.
*           unchanged on exit.
*
*  alpha  - complex         .
*           on entry, alpha specifies the scalar alpha.
*           unchanged on exit.
*
*  a      - complex          array of dimension ( lda, ka ), where ka is
*           k  when  transa = 'n' or 'n',  and is  m  otherwise.
*           before entry with  transa = 'n' or 'n',  the leading  m by k
*           part of the array  a  must contain the matrix  a,  otherwise
*           the leading  k by m  part of the array  a  must contain  the
*           matrix a.
*           unchanged on exit.
*
*  lda    - integer.
*           on entry, lda specifies the first dimension of a as declared
*           in the calling (sub) program. when  transa = 'n' or 'n' then
*           lda must be at least  max( 1, m ), otherwise  lda must be at
*           least  max( 1, k ).
*           unchanged on exit.
*
*  b      - complex          array of dimension ( ldb, kb ), where kb is
*           n  when  transb = 'n' or 'n',  and is  k  otherwise.
*           before entry with  transb = 'n' or 'n',  the leading  k by n
*           part of the array  b  must contain the matrix  b,  otherwise
*           the leading  n by k  part of the array  b  must contain  the
*           matrix b.
*           unchanged on exit.
*
*  ldb    - integer.
*           on entry, ldb specifies the first dimension of b as declared
*           in the calling (sub) program. when  transb = 'n' or 'n' then
*           ldb must be at least  max( 1, k ), otherwise  ldb must be at
*           least  max( 1, n ).
*           unchanged on exit.
*
*  beta   - complex         .
*           on entry,  beta  specifies the scalar  beta.  when  beta  is
*           supplied as zero then c need not be set on input.
*           unchanged on exit.
*
*  c      - complex          array of dimension ( ldc, n ).
*           before entry, the leading  m by n  part of the array  c must
*           contain the matrix  c,  except when  beta  is zero, in which
*           case c need not be set on entry.
*           on exit, the array  c  is overwritten by the  m by n  matrix
*           ( alpha*op( a )*op( b ) + beta*c ).
*
*  ldc    - integer.
*           on entry, ldc specifies the first dimension of c as declared
*           in  the  calling  (sub)  program.   ldc  must  be  at  least
*           max( 1, m ).
*           unchanged on exit.
*
*
*  level 3 blas routine.
*
*  -- written on 8-february-1989.
*     jack dongarra, argonne national laboratory.
*     iain duff, aere harwell.
*     jeremy du croz, numerical algorithms group ltd.
*     sven hammarling, numerical algorithms group ltd.
*
*
*     .. external functions ..
      character*1        dcaptl
      external           dcaptl
*     .. external subroutines ..
      external           lnkerr
*     .. intrinsic functions ..
      intrinsic          conjg, max
*     .. local scalars ..
      logical            conja, conjb, nota, notb
      integer            i, info, j, l, ncola, nrowa, nrowb
      complex *16        temp
*     .. parameters ..
      complex *16        one
      parameter        ( one  = ( 1.0d+0, 0.0d+0 ) )
      complex *16        zero
      parameter        ( zero = ( 0.0d+0, 0.0d+0 ) )
c     .. common blocks ..
      integer inp, iout
      common/io/inp,iout
*     ..
*     .. executable statements ..
*
*     set  nota  and  notb  as  true if  a  and  b  respectively are not
*     conjugated or transposed, set  conja and conjb  as true if  a  and
*     b  respectively are to be  transposed but  not conjugated  and set
*     nrowa, ncola and  nrowb  as the number of rows and  columns  of  a
*     and the number of rows of  b  respectively.
*
      call zgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb,
     $             beta, c, ldc )
*
      return
*
*     end of cgemm .
*
      end


