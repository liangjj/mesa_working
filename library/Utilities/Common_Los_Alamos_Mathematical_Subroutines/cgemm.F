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
      logical lsame
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
#ifdef LOCALBLAS
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
*
*     test the input parameters.
*
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
*
*     quick return if possible.
*
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) )
     $   return
*
*     and when  alpha.eq.zero.
*
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
*
*     start the operations.
*
      if( notb )then
         if( nota )then
*
*           form  c := alpha*a*b + beta*c.
*
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
*
*           form  c := alpha*conjg( a' )*b + beta*c.
*
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
*
*           form  c := alpha*a'*b + beta*c
*
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
*
*           form  c := alpha*a*conjg( b' ) + beta*c.
*
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
*
*           form  c := alpha*a*b'          + beta*c
*
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
*
*           form  c := alpha*conjg( a' )*conjg( b' ) + beta*c.
*
            do 280, j = 1, n
               do 270, i = 1, m
                  temp = zero
                  do 260, l = 1, k
                     temp = temp +
     $                      conjg( a( l, i ) )*conjg( b( j, l ) )
  260             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  270          continue
  280       continue
         else
*
*           form  c := alpha*conjg( a' )*b' + beta*c
*
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
*
*           form  c := alpha*a'*conjg( b' ) + beta*c
*
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
*
*           form  c := alpha*a'*b' + beta*c
*
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
*
      return
*
      end
#endif LOCALBLAS
#ifdef MACHINEBLAS
      call zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,
     $           beta,c,ldc)
      return
      end
#endif MACHINEBLAS


