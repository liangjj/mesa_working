*deck csymm
      subroutine csymm (side, uplo, m, n, alpha, a, lda, b, ldb, beta,
     $   c, ldc)
c***begin prologue  csymm
c***purpose  multiply a complex general matrix by a complex symmetric
c            matrix.
c***library   slatec (blas)
c***category  d1b6
c***type      complex (ssymm-s, dsymm-d, csymm-c)
c***keywords  level 3 blas, linear algebra
c***author  dongarra, j., (anl)
c           duff, i., (aere)
c           du croz, j., (nag)
c           hammarling, s. (nag)
c***description
c
c  csymm  performs one of the matrix-matrix operations
c
c     c := alpha*a*b + beta*c,
c
c  or
c
c     c := alpha*b*a + beta*c,
c
c  where  alpha and beta are scalars, a is a symmetric matrix and  b and
c  c are m by n matrices.
c
c  parameters
c  ==========
c
c  side   - character*1.
c           on entry,  side  specifies whether  the  symmetric matrix  a
c           appears on the  left or right  in the  operation as follows:
c
c              side = 'l' or 'l'   c := alpha*a*b + beta*c,
c
c              side = 'r' or 'r'   c := alpha*b*a + beta*c,
c
c           unchanged on exit.
c
c  uplo   - character*1.
c           on  entry,   uplo  specifies  whether  the  upper  or  lower
c           triangular  part  of  the  symmetric  matrix   a  is  to  be
c           referenced as follows:
c
c              uplo = 'u' or 'u'   only the upper triangular part of the
c                                  symmetric matrix is to be referenced.
c
c              uplo = 'l' or 'l'   only the lower triangular part of the
c                                  symmetric matrix is to be referenced.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry,  m  specifies the number of rows of the matrix  c.
c           m  must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix c.
c           n  must be at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, ka ), where ka is
c           m  when  side = 'l' or 'l'  and is n  otherwise.
c           before entry  with  side = 'l' or 'l',  the  m by m  part of
c           the array  a  must contain the  symmetric matrix,  such that
c           when  uplo = 'u' or 'u', the leading m by m upper triangular
c           part of the array  a  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  a  is not referenced,  and when  uplo = 'l' or 'l',
c           the leading  m by m  lower triangular part  of the  array  a
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  a  is not
c           referenced.
c           before entry  with  side = 'r' or 'r',  the  n by n  part of
c           the array  a  must contain the  symmetric matrix,  such that
c           when  uplo = 'u' or 'u', the leading n by n upper triangular
c           part of the array  a  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  a  is not referenced,  and when  uplo = 'l' or 'l',
c           the leading  n by n  lower triangular part  of the  array  a
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  a  is not
c           referenced.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the  calling (sub) program. when  side = 'l' or 'l'  then
c           lda must be at least  max( 1, m ), otherwise  lda must be at
c           least max( 1, n ).
c           unchanged on exit.
c
c  b      - complex          array of dimension ( ldb, n ).
c           before entry, the leading  m by n part of the array  b  must
c           contain the matrix b.
c           unchanged on exit.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in  the  calling  (sub)  program.   ldb  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c  beta   - complex         .
c           on entry,  beta  specifies the scalar  beta.  when  beta  is
c           supplied as zero then c need not be set on input.
c           unchanged on exit.
c
c  c      - complex          array of dimension ( ldc, n ).
c           before entry, the leading  m by n  part of the array  c must
c           contain the matrix  c,  except when  beta  is zero, in which
c           case c need not be set on entry.
c           on exit, the array  c  is overwritten by the  m by n updated
c           matrix.
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c***references  dongarra, j., du croz, j., duff, i., and hammarling, s.
c                 a set of level 3 basic linear algebra subprograms.
c                 acm toms, vol. 16, no. 1, pp. 1-17, march 1990.
c***routines called  lsame, xerbla
c***revision history  (yymmdd)
c   890208  date written
c   910605  modified to meet slatec prologue standards.  only comment
c           lines were modified.  (bks)
c***end prologue  csymm
c     .. scalar arguments ..
      character*1        side, uplo
      integer            m, n, lda, ldb, ldc
      complex            alpha, beta
c     .. array arguments ..
      complex            a( lda, * ), b( ldb, * ), c( ldc, * )
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     .. local scalars ..
      logical            upper
      integer            i, info, j, k, nrowa
      complex            temp1, temp2
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c***first executable statement  csymm
c
c     set nrowa as the number of rows of a.
c
      if( lsame( side, 'l' ) )then
         nrowa = m
      else
         nrowa = n
      end if
      upper = lsame( uplo, 'u' )
c
c     test the input parameters.
c
      info = 0
      if(      ( .not.lsame( side, 'l' ) ).and.
     $         ( .not.lsame( side, 'r' ) )      )then
         info = 1
      else if( ( .not.upper              ).and.
     $         ( .not.lsame( uplo, 'l' ) )      )then
         info = 2
      else if( m  .lt.0               )then
         info = 3
      else if( n  .lt.0               )then
         info = 4
      else if( lda.lt.max( 1, nrowa ) )then
         info = 7
      else if( ldb.lt.max( 1, m     ) )then
         info = 9
      else if( ldc.lt.max( 1, m     ) )then
         info = 12
      end if
      if( info.ne.0 )then
         call xerbla( 'csymm ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         if( beta.eq.zero )then
            do 20, j = 1, n
               do 10, i = 1, m
                  c( i, j ) = zero
   10          continue
   20       continue
         else
            do 40, j = 1, n
               do 30, i = 1, m
                  c( i, j ) = beta*c( i, j )
   30          continue
   40       continue
         end if
         return
      end if
c
c     start the operations.
c
      if( lsame( side, 'l' ) )then
c
c        form  c := alpha*a*b + beta*c.
c
         if( upper )then
            do 70, j = 1, n
               do 60, i = 1, m
                  temp1 = alpha*b( i, j )
                  temp2 = zero
                  do 50, k = 1, i - 1
                     c( k, j ) = c( k, j ) + temp1    *a( k, i )
                     temp2     = temp2     + b( k, j )*a( k, i )
   50             continue
                  if( beta.eq.zero )then
                     c( i, j ) = temp1*a( i, i ) + alpha*temp2
                  else
                     c( i, j ) = beta *c( i, j ) +
     $                           temp1*a( i, i ) + alpha*temp2
                  end if
   60          continue
   70       continue
         else
            do 100, j = 1, n
               do 90, i = m, 1, -1
                  temp1 = alpha*b( i, j )
                  temp2 = zero
                  do 80, k = i + 1, m
                     c( k, j ) = c( k, j ) + temp1    *a( k, i )
                     temp2     = temp2     + b( k, j )*a( k, i )
   80             continue
                  if( beta.eq.zero )then
                     c( i, j ) = temp1*a( i, i ) + alpha*temp2
                  else
                     c( i, j ) = beta *c( i, j ) +
     $                           temp1*a( i, i ) + alpha*temp2
                  end if
   90          continue
  100       continue
         end if
      else
c
c        form  c := alpha*b*a + beta*c.
c
         do 170, j = 1, n
            temp1 = alpha*a( j, j )
            if( beta.eq.zero )then
               do 110, i = 1, m
                  c( i, j ) = temp1*b( i, j )
  110          continue
            else
               do 120, i = 1, m
                  c( i, j ) = beta*c( i, j ) + temp1*b( i, j )
  120          continue
            end if
            do 140, k = 1, j - 1
               if( upper )then
                  temp1 = alpha*a( k, j )
               else
                  temp1 = alpha*a( j, k )
               end if
               do 130, i = 1, m
                  c( i, j ) = c( i, j ) + temp1*b( i, k )
  130          continue
  140       continue
            do 160, k = j + 1, n
               if( upper )then
                  temp1 = alpha*a( j, k )
               else
                  temp1 = alpha*a( k, j )
               end if
               do 150, i = 1, m
                  c( i, j ) = c( i, j ) + temp1*b( i, k )
  150          continue
  160       continue
  170    continue
      end if
c
      return
c
c     end of csymm .
c
      end
