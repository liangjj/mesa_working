*deck cher
      subroutine cher (uplo, n, alpha, x, incx, a, lda)
c***begin prologue  cher
c***purpose  perform hermitian rank 1 update of a complex hermitian
c            matrix.
c***library   slatec (blas)
c***category  d1b4
c***type      complex (sher-s, dher-d, cher-c)
c***keywords  level 2 blas, linear algebra
c***author  dongarra, j. j., (anl)
c           du croz, j., (nag)
c           hammarling, s., (nag)
c           hanson, r. j., (snla)
c***description
c
c  cher   performs the hermitian rank 1 operation
c
c     a := alpha*x*conjg( x') + a,
c
c  where alpha is a real scalar, x is an n element vector and a is an
c  n by n hermitian matrix.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the array a is to be referenced as
c           follows:
c
c              uplo = 'u' or 'u'   only the upper triangular part of a
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the lower triangular part of a
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular part of the hermitian matrix and the strictly
c           lower triangular part of a is not referenced. on exit, the
c           upper triangular part of the array a is overwritten by the
c           upper triangular part of the updated matrix.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular part of the hermitian matrix and the strictly
c           upper triangular part of a is not referenced. on exit, the
c           lower triangular part of the array a is overwritten by the
c           lower triangular part of the updated matrix.
c           note that the imaginary parts of the diagonal elements need
c           not be set, they are assumed to be zero, and on exit they
c           are set to zero.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c***references  dongarra, j. j., du croz, j., hammarling, s., and
c                 hanson, r. j.  an extended set of fortran basic linear
c                 algebra subprograms.  acm toms, vol. 14, no. 1,
c                 pp. 1-17, march 1988.
c***routines called  lsame, xerbla
c***revision history  (yymmdd)
c   861022  date written
c   910605  modified to meet slatec prologue standards.  only comment
c           lines were modified.  (bks)
c***end prologue  cher
c     .. scalar arguments ..
      real               alpha
      integer            incx, lda, n
      character*1        uplo
c     .. array arguments ..
      complex            a( lda, * ), x( * )
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jx, kx
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max, real
c***first executable statement  cher
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( lda.lt.max( 1, n ) )then
         info = 7
      end if
      if( info.ne.0 )then
         call xerbla( 'cher  ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( alpha.eq.real( zero ) ) )
     $   return
c
c     set the start point in x if the increment is not unity.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through the triangular part
c     of a.
c
      if( lsame( uplo, 'u' ) )then
c
c        form  a  when a is stored in upper triangle.
c
         if( incx.eq.1 )then
            do 20, j = 1, n
               if( x( j ).ne.zero )then
                  temp = alpha*conjg( x( j ) )
                  do 10, i = 1, j - 1
                     a( i, j ) = a( i, j ) + x( i )*temp
   10             continue
                  a( j, j ) = real( a( j, j ) ) + real( x( j )*temp )
               else
                  a( j, j ) = real( a( j, j ) )
               end if
   20       continue
         else
            jx = kx
            do 40, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*conjg( x( jx ) )
                  ix   = kx
                  do 30, i = 1, j - 1
                     a( i, j ) = a( i, j ) + x( ix )*temp
                     ix        = ix        + incx
   30             continue
                  a( j, j ) = real( a( j, j ) ) + real( x( jx )*temp )
               else
                  a( j, j ) = real( a( j, j ) )
               end if
               jx = jx + incx
   40       continue
         end if
      else
c
c        form  a  when a is stored in lower triangle.
c
         if( incx.eq.1 )then
            do 60, j = 1, n
               if( x( j ).ne.zero )then
                  temp      = alpha*conjg( x( j ) )
                  a( j, j ) = real( a( j, j ) ) + real( temp*x( j ) )
                  do 50, i = j + 1, n
                     a( i, j ) = a( i, j ) + x( i )*temp
   50             continue
               else
                  a( j, j ) = real( a( j, j ) )
               end if
   60       continue
         else
            jx = kx
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp      = alpha*conjg( x( jx ) )
                  a( j, j ) = real( a( j, j ) ) + real( temp*x( jx ) )
                  ix        = jx
                  do 70, i = j + 1, n
                     ix        = ix        + incx
                     a( i, j ) = a( i, j ) + x( ix )*temp
   70             continue
               else
                  a( j, j ) = real( a( j, j ) )
               end if
               jx = jx + incx
   80       continue
         end if
      end if
c
      return
c
c     end of cher  .
c
      end
