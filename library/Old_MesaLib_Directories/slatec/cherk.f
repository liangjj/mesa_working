*deck cherk
      subroutine cherk (uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
c***begin prologue  cherk
c***purpose  perform hermitian rank k update of a complex hermitian
c            matrix.
c***library   slatec (blas)
c***category  d1b6
c***type      complex (sherk-s, dherk-d, cherk-c)
c***keywords  level 3 blas, linear algebra
c***author  dongarra, j., (anl)
c           duff, i., (aere)
c           du croz, j., (nag)
c           hammarling, s. (nag)
c***description
c
c  cherk  performs one of the hermitian rank k operations
c
c     c := alpha*a*conjg( a' ) + beta*c,
c
c  or
c
c     c := alpha*conjg( a' )*a + beta*c,
c
c  where  alpha and beta  are  real scalars,  c is an  n by n  hermitian
c  matrix and  a  is an  n by k  matrix in the  first case and a  k by n
c  matrix in the second case.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on  entry,   uplo  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  c  is to be  referenced  as
c           follows:
c
c              uplo = 'u' or 'u'   only the  upper triangular part of  c
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the  lower triangular part of  c
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry,  trans  specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   c := alpha*a*conjg( a' ) + beta*c.
c
c              trans = 'c' or 'c'   c := alpha*conjg( a' )*a + beta*c.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry,  n specifies the order of the matrix c.  n must be
c           at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry with  trans = 'n' or 'n',  k  specifies  the number
c           of  columns   of  the   matrix   a,   and  on   entry   with
c           trans = 'c' or 'c',  k  specifies  the number of rows of the
c           matrix a.  k must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, ka ), where ka is
c           k  when  trans = 'n' or 'n',  and is  n  otherwise.
c           before entry with  trans = 'n' or 'n',  the  leading  n by k
c           part of the array  a  must contain the matrix  a,  otherwise
c           the leading  k by n  part of the array  a  must contain  the
c           matrix a.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in  the  calling  (sub)  program.   when  trans = 'n' or 'n'
c           then  lda must be at least  max( 1, n ), otherwise  lda must
c           be at least  max( 1, k ).
c           unchanged on exit.
c
c  beta   - real            .
c           on entry, beta specifies the scalar beta.
c           unchanged on exit.
c
c  c      - complex          array of dimension ( ldc, n ).
c           before entry  with  uplo = 'u' or 'u',  the leading  n by n
c           upper triangular part of the array c must contain the upper
c           triangular part  of the  hermitian matrix  and the strictly
c           lower triangular part of c is not referenced.  on exit, the
c           upper triangular part of the array  c is overwritten by the
c           upper triangular part of the updated matrix.
c           before entry  with  uplo = 'l' or 'l',  the leading  n by n
c           lower triangular part of the array c must contain the lower
c           triangular part  of the  hermitian matrix  and the strictly
c           upper triangular part of c is not referenced.  on exit, the
c           lower triangular part of the array  c is overwritten by the
c           lower triangular part of the updated matrix.
c           note that the imaginary parts of the diagonal elements need
c           not be set,  they are assumed to be zero,  and on exit they
c           are set to zero.
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, n ).
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
c***end prologue  cherk
c     .. scalar arguments ..
      character*1        uplo, trans
      integer            n, k, lda, ldc
      real               alpha, beta
c     .. array arguments ..
      complex            a( lda, * ), c( ldc, * )
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          cmplx, conjg, max, real
c     .. local scalars ..
      logical            upper
      integer            i, info, j, l, nrowa
      real               rtemp
      complex            temp
c     .. parameters ..
      real               one ,         zero
      parameter        ( one = 1.0e+0, zero = 0.0e+0 )
c***first executable statement  cherk
c
c     test the input parameters.
c
      if( lsame( trans, 'n' ) )then
         nrowa = n
      else
         nrowa = k
      end if
      upper = lsame( uplo, 'u' )
c
      info = 0
      if(      ( .not.upper               ).and.
     $         ( .not.lsame( uplo , 'l' ) )      )then
         info = 1
      else if( ( .not.lsame( trans, 'n' ) ).and.
     $         ( .not.lsame( trans, 'c' ) )      )then
         info = 2
      else if( n  .lt.0               )then
         info = 3
      else if( k  .lt.0               )then
         info = 4
      else if( lda.lt.max( 1, nrowa ) )then
         info = 7
      else if( ldc.lt.max( 1, n     ) )then
         info = 10
      end if
      if( info.ne.0 )then
         call xerbla( 'cherk ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.
     $    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         if( upper )then
            if( beta.eq.zero )then
               do 20, j = 1, n
                  do 10, i = 1, j
                     c( i, j ) = zero
   10             continue
   20          continue
            else
               do 40, j = 1, n
                  do 30, i = 1, j - 1
                     c( i, j ) = beta*c( i, j )
   30             continue
                  c( j, j ) = beta*real( c( j, j ) )
   40          continue
            end if
         else
            if( beta.eq.zero )then
               do 60, j = 1, n
                  do 50, i = j, n
                     c( i, j ) = zero
   50             continue
   60          continue
            else
               do 80, j = 1, n
                  c( j, j ) = beta*real( c( j, j ) )
                  do 70, i = j + 1, n
                     c( i, j ) = beta*c( i, j )
   70             continue
   80          continue
            end if
         end if
         return
      end if
c
c     start the operations.
c
      if( lsame( trans, 'n' ) )then
c
c        form  c := alpha*a*conjg( a' ) + beta*c.
c
         if( upper )then
            do 130, j = 1, n
               if( beta.eq.zero )then
                  do 90, i = 1, j
                     c( i, j ) = zero
   90             continue
               else if( beta.ne.one )then
                  do 100, i = 1, j - 1
                     c( i, j ) = beta*c( i, j )
  100             continue
                  c( j, j ) = beta*real( c( j, j ) )
               end if
               do 120, l = 1, k
                  if( a( j, l ).ne.cmplx( zero ) )then
                     temp = alpha*conjg( a( j, l ) )
                     do 110, i = 1, j - 1
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  110                continue
                     c( j, j ) = real( c( j, j )      ) +
     $                           real( temp*a( i, l ) )
                  end if
  120          continue
  130       continue
         else
            do 180, j = 1, n
               if( beta.eq.zero )then
                  do 140, i = j, n
                     c( i, j ) = zero
  140             continue
               else if( beta.ne.one )then
                  c( j, j ) = beta*real( c( j, j ) )
                  do 150, i = j + 1, n
                     c( i, j ) = beta*c( i, j )
  150             continue
               end if
               do 170, l = 1, k
                  if( a( j, l ).ne.cmplx( zero ) )then
                     temp      = alpha*conjg( a( j, l ) )
                     c( j, j ) = real( c( j, j )      )   +
     $                           real( temp*a( j, l ) )
                     do 160, i = j + 1, n
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  160                continue
                  end if
  170          continue
  180       continue
         end if
      else
c
c        form  c := alpha*conjg( a' )*a + beta*c.
c
         if( upper )then
            do 220, j = 1, n
               do 200, i = 1, j - 1
                  temp = zero
                  do 190, l = 1, k
                     temp = temp + conjg( a( l, i ) )*a( l, j )
  190             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  200          continue
               rtemp = zero
               do 210, l = 1, k
                  rtemp = rtemp + conjg( a( l, j ) )*a( l, j )
  210          continue
               if( beta.eq.zero )then
                  c( j, j ) = alpha*rtemp
               else
                  c( j, j ) = alpha*rtemp + beta*real( c( j, j ) )
               end if
  220       continue
         else
            do 260, j = 1, n
               rtemp = zero
               do 230, l = 1, k
                  rtemp = rtemp + conjg( a( l, j ) )*a( l, j )
  230          continue
               if( beta.eq.zero )then
                  c( j, j ) = alpha*rtemp
               else
                  c( j, j ) = alpha*rtemp + beta*real( c( j, j ) )
               end if
               do 250, i = j + 1, n
                  temp = zero
                  do 240, l = 1, k
                     temp = temp + conjg( a( l, i ) )*a( l, j )
  240             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  250          continue
  260       continue
         end if
      end if
c
      return
c
c     end of cherk .
c
      end
