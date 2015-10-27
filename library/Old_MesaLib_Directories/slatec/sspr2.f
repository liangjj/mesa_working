*deck sspr2
      subroutine sspr2 (uplo, n, alpha, x, incx, y, incy, ap)
c***begin prologue  sspr2
c***purpose  perform the symmetric rank 2 operation.
c***library   slatec (blas)
c***category  d1b4
c***type      single precision (sspr2-s, dspr2-d, cspr2-c)
c***keywords  level 2 blas, linear algebra
c***author  dongarra, j. j., (anl)
c           du croz, j., (nag)
c           hammarling, s., (nag)
c           hanson, r. j., (snla)
c***description
c
c  sspr2  performs the symmetric rank 2 operation
c
c     a := alpha*x*y' + alpha*y*x' + a,
c
c  where alpha is a scalar, x and y are n element vectors and a is an
c  n by n symmetric matrix, supplied in packed form.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the matrix a is supplied in the packed
c           array ap as follows:
c
c              uplo = 'u' or 'u'   the upper triangular part of a is
c                                  supplied in ap.
c
c              uplo = 'l' or 'l'   the lower triangular part of a is
c                                  supplied in ap.
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
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1)*abs( incx)).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  y      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y.
c           unchanged on exit.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c  ap     - real             array of dimension at least
c           ( ( n*( n + 1 ) )/2 ).
c           before entry with  uplo = 'u' or 'u', the array ap must
c           contain the upper triangular part of the symmetric matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
c           and a( 2, 2 ) respectively, and so on. on exit, the array
c           ap is overwritten by the upper triangular part of the
c           updated matrix.
c           before entry with uplo = 'l' or 'l', the array ap must
c           contain the lower triangular part of the symmetric matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
c           and a( 3, 1 ) respectively, and so on. on exit, the array
c           ap is overwritten by the lower triangular part of the
c           updated matrix.
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
c***end prologue  sspr2
c     .. scalar arguments ..
      real               alpha
      integer            incx, incy, n
      character*1        uplo
c     .. array arguments ..
      real               ap( * ), x( * ), y( * )
c     .. parameters ..
      real               zero
      parameter        ( zero = 0.0e+0 )
c     .. local scalars ..
      real               temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, k, kk, kx, ky
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c***first executable statement  sspr2
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
      else if( incy.eq.0 )then
         info = 7
      end if
      if( info.ne.0 )then
         call xerbla( 'sspr2 ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
c
c     set up the start points in x and y if the increments are not both
c     unity.
c
      if( ( incx.ne.1 ).or.( incy.ne.1 ) )then
         if( incx.gt.0 )then
            kx = 1
         else
            kx = 1 - ( n - 1 )*incx
         end if
         if( incy.gt.0 )then
            ky = 1
         else
            ky = 1 - ( n - 1 )*incy
         end if
         jx = kx
         jy = ky
      end if
c
c     start the operations. in this version the elements of the array ap
c     are accessed sequentially with one pass through ap.
c
      kk = 1
      if( lsame( uplo, 'u' ) )then
c
c        form  a  when upper triangle is stored in ap.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 20, j = 1, n
               if( ( x( j ).ne.zero ).or.( y( j ).ne.zero ) )then
                  temp1 = alpha*y( j )
                  temp2 = alpha*x( j )
                  k     = kk
                  do 10, i = 1, j
                     ap( k ) = ap( k ) + x( i )*temp1 + y( i )*temp2
                     k       = k       + 1
   10             continue
               end if
               kk = kk + j
   20       continue
         else
            do 40, j = 1, n
               if( ( x( jx ).ne.zero ).or.( y( jy ).ne.zero ) )then
                  temp1 = alpha*y( jy )
                  temp2 = alpha*x( jx )
                  ix    = kx
                  iy    = ky
                  do 30, k = kk, kk + j - 1
                     ap( k ) = ap( k ) + x( ix )*temp1 + y( iy )*temp2
                     ix      = ix      + incx
                     iy      = iy      + incy
   30             continue
               end if
               jx = jx + incx
               jy = jy + incy
               kk = kk + j
   40       continue
         end if
      else
c
c        form  a  when lower triangle is stored in ap.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               if( ( x( j ).ne.zero ).or.( y( j ).ne.zero ) )then
                  temp1 = alpha*y( j )
                  temp2 = alpha*x( j )
                  k     = kk
                  do 50, i = j, n
                     ap( k ) = ap( k ) + x( i )*temp1 + y( i )*temp2
                     k       = k       + 1
   50             continue
               end if
               kk = kk + n - j + 1
   60       continue
         else
            do 80, j = 1, n
               if( ( x( jx ).ne.zero ).or.( y( jy ).ne.zero ) )then
                  temp1 = alpha*y( jy )
                  temp2 = alpha*x( jx )
                  ix    = jx
                  iy    = jy
                  do 70, k = kk, kk + n - j
                     ap( k ) = ap( k ) + x( ix )*temp1 + y( iy )*temp2
                     ix      = ix      + incx
                     iy      = iy      + incy
   70             continue
               end if
               jx = jx + incx
               jy = jy + incy
               kk = kk + n - j + 1
   80       continue
         end if
      end if
c
      return
c
c     end of sspr2 .
c
      end
