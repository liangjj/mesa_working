*deck ssymv
      subroutine ssymv (uplo, n, alpha, a, lda, x, incx, beta, y, incy)
c***begin prologue  ssymv
c***purpose  multiply a real vector by a real symmetric matrix.
c***library   slatec (blas)
c***category  d1b4
c***type      single precision (ssymv-s, dsymv-d, csymv-c)
c***keywords  level 2 blas, linear algebra
c***author  dongarra, j. j., (anl)
c           du croz, j., (nag)
c           hammarling, s., (nag)
c           hanson, r. j., (snla)
c***description
c
c  ssymv  performs the matrix-vector  operation
c
c     y := alpha*a*x + beta*y,
c
c  where alpha and beta are scalars, x and y are n element vectors and
c  a is an n by n symmetric matrix.
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
c  a      - real             array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular part of the symmetric matrix and the strictly
c           lower triangular part of a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular part of the symmetric matrix and the strictly
c           upper triangular part of a is not referenced.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c  x      - real             array of dimension at least
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
c  beta   - real            .
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.
c
c  y      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y. on exit, y is overwritten by the updated
c           vector y.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
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
c***end prologue  ssymv
c     .. scalar arguments ..
      real               alpha, beta
      integer            incx, incy, lda, n
      character*1        uplo
c     .. array arguments ..
      real               a( lda, * ), x( * ), y( * )
c     .. parameters ..
      real               one         , zero
      parameter        ( one = 1.0e+0, zero = 0.0e+0 )
c     .. local scalars ..
      real               temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, kx, ky
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c***first executable statement  ssymv
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( lda.lt.max( 1, n ) )then
         info = 5
      else if( incx.eq.0 )then
         info = 7
      else if( incy.eq.0 )then
         info = 10
      end if
      if( info.ne.0 )then
         call xerbla( 'ssymv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     set up the start points in  x  and  y.
c
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
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through the triangular part
c     of a.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, n
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, n
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, n
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, n
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      if( lsame( uplo, 'u' ) )then
c
c        form  y  when a is stored in upper triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               temp1 = alpha*x( j )
               temp2 = zero
               do 50, i = 1, j - 1
                  y( i ) = y( i ) + temp1*a( i, j )
                  temp2  = temp2  + a( i, j )*x( i )
   50          continue
               y( j ) = y( j ) + temp1*a( j, j ) + alpha*temp2
   60       continue
         else
            jx = kx
            jy = ky
            do 80, j = 1, n
               temp1 = alpha*x( jx )
               temp2 = zero
               ix    = kx
               iy    = ky
               do 70, i = 1, j - 1
                  y( iy ) = y( iy ) + temp1*a( i, j )
                  temp2   = temp2   + a( i, j )*x( ix )
                  ix      = ix      + incx
                  iy      = iy      + incy
   70          continue
               y( jy ) = y( jy ) + temp1*a( j, j ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
   80       continue
         end if
      else
c
c        form  y  when a is stored in lower triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 100, j = 1, n
               temp1  = alpha*x( j )
               temp2  = zero
               y( j ) = y( j )       + temp1*a( j, j )
               do 90, i = j + 1, n
                  y( i ) = y( i ) + temp1*a( i, j )
                  temp2  = temp2  + a( i, j )*x( i )
   90          continue
               y( j ) = y( j ) + alpha*temp2
  100       continue
         else
            jx = kx
            jy = ky
            do 120, j = 1, n
               temp1   = alpha*x( jx )
               temp2   = zero
               y( jy ) = y( jy )       + temp1*a( j, j )
               ix      = jx
               iy      = jy
               do 110, i = j + 1, n
                  ix      = ix      + incx
                  iy      = iy      + incy
                  y( iy ) = y( iy ) + temp1*a( i, j )
                  temp2   = temp2   + a( i, j )*x( ix )
  110          continue
               y( jy ) = y( jy ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
  120       continue
         end if
      end if
c
      return
c
c     end of ssymv .
c
      end
