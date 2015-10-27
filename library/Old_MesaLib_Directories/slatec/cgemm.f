*deck cgemm
      subroutine cgemm (transa, transb, m, n, k, alpha, a, lda, b, ldb,
     $   beta, c, ldc)
c***begin prologue  cgemm
c***purpose  multiply a complex general matrix by a complex general
c            matrix.
c***library   slatec (blas)
c***category  d1b6
c***type      complex (sgemm-s, dgemm-d, cgemm-c)
c***keywords  level 3 blas, linear algebra
c***author  dongarra, j., (anl)
c           duff, i., (aere)
c           du croz, j., (nag)
c           hammarling, s. (nag)
c***description
c
c  cgemm  performs one of the matrix-matrix operations
c
c     c := alpha*op( a )*op( b ) + beta*c,
c
c  where  op( x ) is one of
c
c     op( x ) = x   or   op( x ) = x'   or   op( x ) = conjg( x' ),
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
c              transa = 'c' or 'c',  op( a ) = conjg( a' ).
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
c              transb = 'c' or 'c',  op( b ) = conjg( b' ).
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
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, ka ), where ka is
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
c  b      - complex          array of dimension ( ldb, kb ), where kb is
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
c  beta   - complex         .
c           on entry,  beta  specifies the scalar  beta.  when  beta  is
c           supplied as zero then c need not be set on input.
c           unchanged on exit.
c
c  c      - complex          array of dimension ( ldc, n ).
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
c***references  dongarra, j., du croz, j., duff, i., and hammarling, s.
c                 a set of level 3 basic linear algebra subprograms.
c                 acm toms, vol. 16, no. 1, pp. 1-17, march 1990.
c***routines called  lsame, xerbla
c***revision history  (yymmdd)
c   890208  date written
c   910605  modified to meet slatec prologue standards.  only comment
c           lines were modified.  (bks)
c***end prologue  cgemm
c     .. scalar arguments ..
      character*1        transa, transb
      integer            m, n, k, lda, ldb, ldc
      complex            alpha, beta
c     .. array arguments ..
      complex            a( lda, * ), b( ldb, * ), c( ldc, * )
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max
c     .. local scalars ..
      logical            conja, conjb, nota, notb
      integer            i, info, j, l, ncola, nrowa, nrowb
      complex            temp
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c***first executable statement  cgemm
c
c     set  nota  and  notb  as  true if  a  and  b  respectively are not
c     conjugated or transposed, set  conja and conjb  as true if  a  and
c     b  respectively are to be  transposed but  not conjugated  and set
c     nrowa, ncola and  nrowb  as the number of rows and  columns  of  a
c     and the number of rows of  b  respectively.
c
      nota  = lsame( transa, 'n' )
      notb  = lsame( transb, 'n' )
      conja = lsame( transa, 'c' )
      conjb = lsame( transb, 'c' )
      if( nota )then
         nrowa = m
         ncola = k
      else
         nrowa = k
         ncola = m
      end if
      if( notb )then
         nrowb = k
      else
         nrowb = n
      end if
c
c     test the input parameters.
c
      info = 0
      if(      ( .not.nota                 ).and.
     $         ( .not.conja                ).and.
     $         ( .not.lsame( transa, 't' ) )      )then
         info = 1
      else if( ( .not.notb                 ).and.
     $         ( .not.conjb                ).and.
     $         ( .not.lsame( transb, 't' ) )      )then
         info = 2
      else if( m  .lt.0               )then
         info = 3
      else if( n  .lt.0               )then
         info = 4
      else if( k  .lt.0               )then
         info = 5
      else if( lda.lt.max( 1, nrowa ) )then
         info = 8
      else if( ldb.lt.max( 1, nrowb ) )then
         info = 10
      else if( ldc.lt.max( 1, m     ) )then
         info = 13
      end if
      if( info.ne.0 )then
         call xerbla( 'cgemm ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) )
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
      if( notb )then
         if( nota )then
c
c           form  c := alpha*a*b + beta*c.
c
            do 90, j = 1, n
               if( beta.eq.zero )then
                  do 50, i = 1, m
                     c( i, j ) = zero
   50             continue
               else if( beta.ne.one )then
                  do 60, i = 1, m
                     c( i, j ) = beta*c( i, j )
   60             continue
               end if
               do 80, l = 1, k
                  if( b( l, j ).ne.zero )then
                     temp = alpha*b( l, j )
                     do 70, i = 1, m
                        c( i, j ) = c( i, j ) + temp*a( i, l )
   70                continue
                  end if
   80          continue
   90       continue
         else if( conja )then
c
c           form  c := alpha*conjg( a' )*b + beta*c.
c
            do 120, j = 1, n
               do 110, i = 1, m
                  temp = zero
                  do 100, l = 1, k
                     temp = temp + conjg( a( l, i ) )*b( l, j )
  100             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  110          continue
  120       continue
         else
c
c           form  c := alpha*a'*b + beta*c
c
            do 150, j = 1, n
               do 140, i = 1, m
                  temp = zero
                  do 130, l = 1, k
                     temp = temp + a( l, i )*b( l, j )
  130             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  140          continue
  150       continue
         end if
      else if( nota )then
         if( conjb )then
c
c           form  c := alpha*a*conjg( b' ) + beta*c.
c
            do 200, j = 1, n
               if( beta.eq.zero )then
                  do 160, i = 1, m
                     c( i, j ) = zero
  160             continue
               else if( beta.ne.one )then
                  do 170, i = 1, m
                     c( i, j ) = beta*c( i, j )
  170             continue
               end if
               do 190, l = 1, k
                  if( b( j, l ).ne.zero )then
                     temp = alpha*conjg( b( j, l ) )
                     do 180, i = 1, m
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  180                continue
                  end if
  190          continue
  200       continue
         else
c
c           form  c := alpha*a*b'          + beta*c
c
            do 250, j = 1, n
               if( beta.eq.zero )then
                  do 210, i = 1, m
                     c( i, j ) = zero
  210             continue
               else if( beta.ne.one )then
                  do 220, i = 1, m
                     c( i, j ) = beta*c( i, j )
  220             continue
               end if
               do 240, l = 1, k
                  if( b( j, l ).ne.zero )then
                     temp = alpha*b( j, l )
                     do 230, i = 1, m
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  230                continue
                  end if
  240          continue
  250       continue
         end if
      else if( conja )then
         if( conjb )then
c
c           form  c := alpha*conjg( a' )*conjg( b' ) + beta*c.
c
            do 280, j = 1, n
               do 270, i = 1, m
                  temp = zero
                  do 260, l = 1, k
                     temp = temp + conjg( a( l, i ) )*conjg( b( j, l ) )
  260             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  270          continue
  280       continue
         else
c
c           form  c := alpha*conjg( a' )*b' + beta*c
c
            do 310, j = 1, n
               do 300, i = 1, m
                  temp = zero
                  do 290, l = 1, k
                     temp = temp + conjg( a( l, i ) )*b( j, l )
  290             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  300          continue
  310       continue
         end if
      else
         if( conjb )then
c
c           form  c := alpha*a'*conjg( b' ) + beta*c
c
            do 340, j = 1, n
               do 330, i = 1, m
                  temp = zero
                  do 320, l = 1, k
                     temp = temp + a( l, i )*conjg( b( j, l ) )
  320             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  330          continue
  340       continue
         else
c
c           form  c := alpha*a'*b' + beta*c
c
            do 370, j = 1, n
               do 360, i = 1, m
                  temp = zero
                  do 350, l = 1, k
                     temp = temp + a( l, i )*b( j, l )
  350             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  360          continue
  370       continue
         end if
      end if
c
      return
c
c     end of cgemm .
c
      end
