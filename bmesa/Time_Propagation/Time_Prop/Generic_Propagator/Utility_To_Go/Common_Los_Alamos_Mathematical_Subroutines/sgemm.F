*deck @(#)sgemm.f	5.1 11/29/95
      subroutine sgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb,
     $                   beta, c, ldc )
c***begin prologue     sgemm.f
c***date written       890208  
c***revision date      11/29/95      
c
c***keywords           
c***author          
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c***source             @(#)sgemm.f	5.1   11/29/95
c***purpose            
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
c***description
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
c    
c
c***references
c
c***routines called
c
c***end prologue       sgemm.f
      implicit none
c     .. scalar arguments ..
      character*1        transa, transb
      integer            m, n, k, lda, ldb, ldc
      real*8             alpha, beta
c     .. array arguments ..
      real*8             a( lda, * ), b( ldb, * ), c( ldc, * )
c     ..
c     .. external functions ..
      character*1        dcaptl
      external           dcaptl
c     .. external subroutines ..
      external           lnkerr
c     .. intrinsic functions ..
      intrinsic          max
c     .. local scalars ..
      logical            nota, notb
      logical lsame
      integer            i, info, j, l, ncola, nrowa, nrowb
      real*8             temp
c     .. parameters ..
      real*8             one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     .. common blocks ..
      integer inp,iout
      common/io/inp,iout
c     ..
c     .. executable statements ..
#ifdef LOCALBLAS
*
*     set  nota  and  notb  as  true if  a  and  b  respectively are not
*     transposed and set  nrowa, ncola and  nrowb  as the number of rows
*     and  columns of  a  and the  number of  rows  of  b  respectively.
*
      nota  = lsame( transa, 'n' )
      notb  = lsame( transb, 'n' )
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
     $         ( .not.lsame( transa, 'c' ) ).and.
     $         ( .not.lsame( transa, 't' ) )      )then
         info = 1
      else if( ( .not.notb                 ).and.
     $         ( .not.lsame( transb, 'c' ) ).and.
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
         call xerbla( 'dgemm ', info )
         return
      end if
*
*     quick return if possible.
*
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) )
     $   return
*
*     and if  alpha.eq.zero.
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
         else
*
*           form  c := alpha*a'*b + beta*c
*
            do 120, j = 1, n
               do 110, i = 1, m
                  temp = zero
                  do 100, l = 1, k
                     temp = temp + a( l, i )*b( l, j )
  100             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  110          continue
  120       continue
         end if
      else
         if( nota )then
*
*           form  c := alpha*a*b' + beta*c
*
            do 170, j = 1, n
               if( beta.eq.zero )then
                  do 130, i = 1, m
                     c( i, j ) = zero
  130             continue
               else if( beta.ne.one )then
                  do 140, i = 1, m
                     c( i, j ) = beta*c( i, j )
  140             continue
               end if
               do 160, l = 1, k
                  if( b( j, l ).ne.zero )then
                     temp = alpha*b( j, l )
                     do 150, i = 1, m
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  150                continue
                  end if
  160          continue
  170       continue
         else
*
*           form  c := alpha*a'*b' + beta*c
*
            do 200, j = 1, n
               do 190, i = 1, m
                  temp = zero
                  do 180, l = 1, k
                     temp = temp + a( l, i )*b( j, l )
  180             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  190          continue
  200       continue
         end if
      end if
      return
      end

#endif LOCALBLAS
#ifdef MACHINEBLAS
      call dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
      return
      end
#endif MACHINEBLAS
