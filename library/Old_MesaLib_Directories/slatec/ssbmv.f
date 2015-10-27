*deck ssbmv
      subroutine ssbmv (uplo, n, k, alpha, a, lda, x, incx, beta, y,
     $   incy)
c***begin prologue  ssbmv
c***purpose  multiply a real vector by a real symmetric band matrix.
c***library   slatec (blas)
c***category  d1b4
c***type      single precision (ssbmv-s, dsbmv-d, csbmv-c)
c***keywords  level 2 blas, linear algebra
c***author  dongarra, j. j., (anl)
c           du croz, j., (nag)
c           hammarling, s., (nag)
c           hanson, r. j., (snla)
c***description
c
c  ssbmv  performs the matrix-vector  operation
c
c     y := alpha*a*x + beta*y,
c
c  where alpha and beta are scalars, x and y are n element vectors and
c  a is an n by n symmetric band matrix, with k super-diagonals.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the band matrix a is being supplied as
c           follows:
c
c              uplo = 'u' or 'u'   the upper triangular part of a is
c                                  being supplied.
c
c              uplo = 'l' or 'l'   the lower triangular part of a is
c                                  being supplied.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry, k specifies the number of super-diagonals of the
c           matrix a. k must satisfy  0 .le. k.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry with uplo = 'u' or 'u', the leading ( k + 1 )
c           by n part of the array a must contain the upper triangular
c           band part of the symmetric matrix, supplied column by
c           column, with the leading diagonal of the matrix in row
c           ( k + 1 ) of the array, the first super-diagonal starting at
c           position 2 in row k, and so on. the top left k by k triangle
c           of the array a is not referenced.
c           the following program segment will transfer the upper
c           triangular part of a symmetric band matrix from conventional
c           full matrix storage to band storage:
c
c                 do 20, j = 1, n
c                    m = k + 1 - j
c                    do 10, i = max( 1, j - k ), j
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           before entry with uplo = 'l' or 'l', the leading ( k + 1 )
c           by n part of the array a must contain the lower triangular
c           band part of the symmetric matrix, supplied column by
c           column, with the leading diagonal of the matrix in row 1 of
c           the array, the first sub-diagonal starting at position 1 in
c           row 2, and so on. the bottom right k by k triangle of the
c           array a is not referenced.
c           the following program segment will transfer the lower
c           triangular part of a symmetric band matrix from conventional
c           full matrix storage to band storage:
c
c                 do 20, j = 1, n
c                    m = 1 - j
c                    do 10, i = j, min( n, j + k )
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           ( k + 1 ).
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the
c           vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - real            .
c           on entry, beta specifies the scalar beta.
c           unchanged on exit.
c
c  y      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the
c           vector y. on exit, y is overwritten by the updated vector y.
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
c***end prologue  ssbmv
c     .. scalar arguments ..
      real               alpha, beta
      integer            incx, incy, k, lda, n
      character*1        uplo
c     .. array arguments ..
      real               a( lda, * ), x( * ), y( * )
c     .. parameters ..
      real               one         , zero
      parameter        ( one = 1.0e+0, zero = 0.0e+0 )
c     .. local scalars ..
      real               temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max, min
c***first executable statement  ssbmv
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( k.lt.0 )then
         info = 3
      else if( lda.lt.( k + 1 ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      else if( incy.eq.0 )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'ssbmv ', info )
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
c     start the operations. in this version the elements of the array a
c     are accessed sequentially with one pass through a.
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
c        form  y  when upper triangle of a is stored.
c
         kplus1 = k + 1
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               temp1 = alpha*x( j )
               temp2 = zero
               l     = kplus1 - j
               do 50, i = max( 1, j - k ), j - 1
                  y( i ) = y( i ) + temp1*a( l + i, j )
                  temp2  = temp2  + a( l + i, j )*x( i )
   50          continue
               y( j ) = y( j ) + temp1*a( kplus1, j ) + alpha*temp2
   60       continue
         else
            jx = kx
            jy = ky
            do 80, j = 1, n
               temp1 = alpha*x( jx )
               temp2 = zero
               ix    = kx
               iy    = ky
               l     = kplus1 - j
               do 70, i = max( 1, j - k ), j - 1
                  y( iy ) = y( iy ) + temp1*a( l + i, j )
                  temp2   = temp2   + a( l + i, j )*x( ix )
                  ix      = ix      + incx
                  iy      = iy      + incy
   70          continue
               y( jy ) = y( jy ) + temp1*a( kplus1, j ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
               if( j.gt.k )then
                  kx = kx + incx
                  ky = ky + incy
               end if
   80       continue
         end if
      else
c
c        form  y  when lower triangle of a is stored.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 100, j = 1, n
               temp1  = alpha*x( j )
               temp2  = zero
               y( j ) = y( j )       + temp1*a( 1, j )
               l      = 1            - j
               do 90, i = j + 1, min( n, j + k )
                  y( i ) = y( i ) + temp1*a( l + i, j )
                  temp2  = temp2  + a( l + i, j )*x( i )
   90          continue
               y( j ) = y( j ) + alpha*temp2
  100       continue
         else
            jx = kx
            jy = ky
            do 120, j = 1, n
               temp1   = alpha*x( jx )
               temp2   = zero
               y( jy ) = y( jy )       + temp1*a( 1, j )
               l       = 1             - j
               ix      = jx
               iy      = jy
               do 110, i = j + 1, min( n, j + k )
                  ix      = ix      + incx
                  iy      = iy      + incy
                  y( iy ) = y( iy ) + temp1*a( l + i, j )
                  temp2   = temp2   + a( l + i, j )*x( ix )
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
c     end of ssbmv .
c
      end
