*deck cgbmv
      subroutine cgbmv (trans, m, n, kl, ku, alpha, a, lda, x, incx,
     $   beta, y, incy)
c***begin prologue  cgbmv
c***purpose  multiply a complex vector by a complex general band matrix.
c***library   slatec (blas)
c***category  d1b4
c***type      complex (sgbmv-s, dgbmv-d, cgbmv-c)
c***keywords  level 2 blas, linear algebra
c***author  dongarra, j. j., (anl)
c           du croz, j., (nag)
c           hammarling, s., (nag)
c           hanson, r. j., (snla)
c***description
c
c  cgbmv  performs one of the matrix-vector operations
c
c     y := alpha*a*x + beta*y,   or   y := alpha*a'*x + beta*y,   or
c
c     y := alpha*conjg( a' )*x + beta*y,
c
c  where alpha and beta are scalars, x and y are vectors and a is an
c  m by n band matrix, with kl sub-diagonals and ku super-diagonals.
c
c  parameters
c  ==========
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   y := alpha*a*x + beta*y.
c
c              trans = 't' or 't'   y := alpha*a'*x + beta*y.
c
c              trans = 'c' or 'c'   y := alpha*conjg( a' )*x + beta*y.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of the matrix a.
c           m must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  kl     - integer.
c           on entry, kl specifies the number of sub-diagonals of the
c           matrix a. kl must satisfy  0 .le. kl.
c           unchanged on exit.
c
c  ku     - integer.
c           on entry, ku specifies the number of super-diagonals of the
c           matrix a. ku must satisfy  0 .le. ku.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry, the leading ( kl + ku + 1 ) by n part of the
c           array a must contain the matrix of coefficients, supplied
c           column by column, with the leading diagonal of the matrix in
c           row ( ku + 1 ) of the array, the first super-diagonal
c           starting at position 2 in row ku, the first sub-diagonal
c           starting at position 1 in row ( ku + 2 ), and so on.
c           elements in the array a that do not correspond to elements
c           in the band matrix (such as the top left ku by ku triangle)
c           are not referenced.
c           the following program segment will transfer a band matrix
c           from conventional full matrix storage to band storage:
c
c                 do 20, j = 1, n
c                    k = ku + 1 - j
c                    do 10, i = max( 1, j - ku ), min( m, j + kl )
c                       a( k + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           ( kl + ku + 1 ).
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
c           before entry, the incremented array x must contain the
c           vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - complex         .
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.
c
c  y      - complex          array of dimension at least
c           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
c           before entry, the incremented array y must contain the
c           vector y. on exit, y is overwritten by the updated vector y.
c
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
c***end prologue  cgbmv
c     .. scalar arguments ..
      complex            alpha, beta
      integer            incx, incy, kl, ku, lda, m, n
      character*1        trans
c     .. array arguments ..
      complex            a( lda, * ), x( * ), y( * )
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, iy, j, jx, jy, k, kup1, kx, ky,
     $                   lenx, leny
      logical            noconj
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max, min
c***first executable statement  cgbmv
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 1
      else if( m.lt.0 )then
         info = 2
      else if( n.lt.0 )then
         info = 3
      else if( kl.lt.0 )then
         info = 4
      else if( ku.lt.0 )then
         info = 5
      else if( lda.lt.( kl + ku + 1 ) )then
         info = 8
      else if( incx.eq.0 )then
         info = 10
      else if( incy.eq.0 )then
         info = 13
      end if
      if( info.ne.0 )then
         call xerbla( 'cgbmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
      noconj = lsame( trans, 't' )
c
c     set  lenx  and  leny, the lengths of the vectors x and y, and set
c     up the start points in  x  and  y.
c
      if( lsame( trans, 'n' ) )then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( lenx - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( leny - 1 )*incy
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through the band part of a.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, leny
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, leny
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, leny
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, leny
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      kup1 = ku + 1
      if( lsame( trans, 'n' ) )then
c
c        form  y := alpha*a*x + y.
c
         jx = kx
         if( incy.eq.1 )then
            do 60, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  k    = kup1 - j
                  do 50, i = max( 1, j - ku ), min( m, j + kl )
                     y( i ) = y( i ) + temp*a( k + i, j )
   50             continue
               end if
               jx = jx + incx
   60       continue
         else
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  iy   = ky
                  k    = kup1 - j
                  do 70, i = max( 1, j - ku ), min( m, j + kl )
                     y( iy ) = y( iy ) + temp*a( k + i, j )
                     iy      = iy      + incy
   70             continue
               end if
               jx = jx + incx
               if( j.gt.ku )
     $            ky = ky + incy
   80       continue
         end if
      else
c
c        form  y := alpha*a'*x + y  or  y := alpha*conjg( a' )*x + y.
c
         jy = ky
         if( incx.eq.1 )then
            do 110, j = 1, n
               temp = zero
               k    = kup1 - j
               if( noconj )then
                  do 90, i = max( 1, j - ku ), min( m, j + kl )
                     temp = temp + a( k + i, j )*x( i )
   90             continue
               else
                  do 100, i = max( 1, j - ku ), min( m, j + kl )
                     temp = temp + conjg( a( k + i, j ) )*x( i )
  100             continue
               end if
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
  110       continue
         else
            do 140, j = 1, n
               temp = zero
               ix   = kx
               k    = kup1 - j
               if( noconj )then
                  do 120, i = max( 1, j - ku ), min( m, j + kl )
                     temp = temp + a( k + i, j )*x( ix )
                     ix   = ix   + incx
  120             continue
               else
                  do 130, i = max( 1, j - ku ), min( m, j + kl )
                     temp = temp + conjg( a( k + i, j ) )*x( ix )
                     ix   = ix   + incx
  130             continue
               end if
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
               if( j.gt.ku )
     $            kx = kx + incx
  140       continue
         end if
      end if
c
      return
c
c     end of cgbmv .
c
      end
