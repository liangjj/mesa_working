*deck sspmv
      subroutine sspmv (uplo, n, alpha, ap, x, incx, beta, y, incy)
c***begin prologue  sspmv
c***purpose  perform the matrix-vector operation.
c***library   slatec (blas)
c***category  d1b4
c***type      single precision (sspmv-s, dspmv-d, cspmv-c)
c***keywords  level 2 blas, linear algebra
c***author  dongarra, j. j., (anl)
c           du croz, j., (nag)
c           hammarling, s., (nag)
c           hanson, r. j., (snla)
c***description
c
c  sspmv  performs the matrix-vector operation
c
c     y := alpha*a*x + beta*y,
c
c  where alpha and beta are scalars, x and y are n element vectors and
c  a is an n by n symmetric matrix, supplied in packed form.
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
c  ap     - real             array of dimension at least
c           ( ( n*( n + 1))/2).
c           before entry with uplo = 'u' or 'u', the array ap must
c           contain the upper triangular part of the symmetric matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
c           and a( 2, 2 ) respectively, and so on.
c           before entry with uplo = 'l' or 'l', the array ap must
c           contain the lower triangular part of the symmetric matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
c           and a( 3, 1 ) respectively, and so on.
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
c***end prologue  sspmv
c     .. scalar arguments ..
      real               alpha, beta
      integer            incx, incy, n
      character*1        uplo
c     .. array arguments ..
      real               ap( * ), x( * ), y( * )
c     .. parameters ..
      real               one         , zero
      parameter        ( one = 1.0e+0, zero = 0.0e+0 )
c     .. local scalars ..
      real               temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, k, kk, kx, ky
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c***first executable statement  sspmv
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
         info = 6
      else if( incy.eq.0 )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'sspmv ', info )
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
c     start the operations. in this version the elements of the array ap
c     are accessed sequentially with one pass through ap.
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
      kk = 1
      if( lsame( uplo, 'u' ) )then
c
c        form  y  when ap contains the upper triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               temp1 = alpha*x( j )
               temp2 = zero
               k     = kk
               do 50, i = 1, j - 1
                  y( i ) = y( i ) + temp1*ap( k )
                  temp2  = temp2  + ap( k )*x( i )
                  k      = k      + 1
   50          continue
               y( j ) = y( j ) + temp1*ap( kk + j - 1 ) + alpha*temp2
               kk     = kk     + j
   60       continue
         else
            jx = kx
            jy = ky
            do 80, j = 1, n
               temp1 = alpha*x( jx )
               temp2 = zero
               ix    = kx
               iy    = ky
               do 70, k = kk, kk + j - 2
                  y( iy ) = y( iy ) + temp1*ap( k )
                  temp2   = temp2   + ap( k )*x( ix )
                  ix      = ix      + incx
                  iy      = iy      + incy
   70          continue
               y( jy ) = y( jy ) + temp1*ap( kk + j - 1 ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
               kk      = kk      + j
   80       continue
         end if
      else
c
c        form  y  when ap contains the lower triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 100, j = 1, n
               temp1  = alpha*x( j )
               temp2  = zero
               y( j ) = y( j )       + temp1*ap( kk )
               k      = kk           + 1
               do 90, i = j + 1, n
                  y( i ) = y( i ) + temp1*ap( k )
                  temp2  = temp2  + ap( k )*x( i )
                  k      = k      + 1
   90          continue
               y( j ) = y( j ) + alpha*temp2
               kk     = kk     + ( n - j + 1 )
  100       continue
         else
            jx = kx
            jy = ky
            do 120, j = 1, n
               temp1   = alpha*x( jx )
               temp2   = zero
               y( jy ) = y( jy )       + temp1*ap( kk )
               ix      = jx
               iy      = jy
               do 110, k = kk + 1, kk + n - j
                  ix      = ix      + incx
                  iy      = iy      + incy
                  y( iy ) = y( iy ) + temp1*ap( k )
                  temp2   = temp2   + ap( k )*x( ix )
  110          continue
               y( jy ) = y( jy ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
               kk      = kk      + ( n - j + 1 )
  120       continue
         end if
      end if
c
      return
c
c     end of sspmv .
c
      end
