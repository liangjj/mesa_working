*deck ctpmv
      subroutine ctpmv (uplo, trans, diag, n, ap, x, incx)
c***begin prologue  ctpmv
c***purpose  perform one of the matrix-vector operations.
c***library   slatec (blas)
c***category  d1b4
c***type      complex (stpmv-s, dtpmv-d, ctpmv-c)
c***keywords  level 2 blas, linear algebra
c***author  dongarra, j. j., (anl)
c           du croz, j., (nag)
c           hammarling, s., (nag)
c           hanson, r. j., (snla)
c***description
c
c  ctpmv  performs one of the matrix-vector operations
c
c     x := a*x,   or   x := a'*x,   or   x := conjg( a')*x,
c
c  where x is an n element vector and  a is an n by n unit, or non-unit,
c  upper or lower triangular matrix, supplied in packed form.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   x := a*x.
c
c              trans = 't' or 't'   x := a'*x.
c
c              trans = 'c' or 'c'   x := conjg( a' )*x.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  ap     - complex          array of dimension at least
c           ( ( n*( n + 1 ) )/2 ).
c           before entry with  uplo = 'u' or 'u', the array ap must
c           contain the upper triangular matrix packed sequentially,
c           column by column, so that ap( 1 ) contains a( 1, 1 ),
c           ap( 2 ) and ap( 3 ) contain a( 1, 2 ) and a( 2, 2 )
c           respectively, and so on.
c           before entry with uplo = 'l' or 'l', the array ap must
c           contain the lower triangular matrix packed sequentially,
c           column by column, so that ap( 1 ) contains a( 1, 1 ),
c           ap( 2 ) and ap( 3 ) contain a( 2, 1 ) and a( 3, 1 )
c           respectively, and so on.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced, but are assumed to be unity.
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x. on exit, x is overwritten with the
c           transformed vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
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
c***end prologue  ctpmv
c     .. scalar arguments ..
      integer            incx, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      complex            ap( * ), x( * )
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jx, k, kk, kx
      logical            noconj, nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg
c***first executable statement  ctpmv
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( incx.eq.0 )then
         info = 7
      end if
      if( info.ne.0 )then
         call xerbla( 'ctpmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      noconj = lsame( trans, 't' )
      nounit = lsame( diag , 'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of ap are
c     accessed sequentially with one pass through ap.
c
      if( lsame( trans, 'n' ) )then
c
c        form  x:= a*x.
c
         if( lsame( uplo, 'u' ) )then
            kk = 1
            if( incx.eq.1 )then
               do 20, j = 1, n
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     k    = kk
                     do 10, i = 1, j - 1
                        x( i ) = x( i ) + temp*ap( k )
                        k      = k      + 1
   10                continue
                     if( nounit )
     $                  x( j ) = x( j )*ap( kk + j - 1 )
                  end if
                  kk = kk + j
   20          continue
            else
               jx = kx
               do 40, j = 1, n
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 30, k = kk, kk + j - 2
                        x( ix ) = x( ix ) + temp*ap( k )
                        ix      = ix      + incx
   30                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*ap( kk + j - 1 )
                  end if
                  jx = jx + incx
                  kk = kk + j
   40          continue
            end if
         else
            kk = ( n*( n + 1 ) )/2
            if( incx.eq.1 )then
               do 60, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     k    = kk
                     do 50, i = n, j + 1, -1
                        x( i ) = x( i ) + temp*ap( k )
                        k      = k      - 1
   50                continue
                     if( nounit )
     $                  x( j ) = x( j )*ap( kk - n + j )
                  end if
                  kk = kk - ( n - j + 1 )
   60          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 80, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 70, k = kk, kk - ( n - ( j + 1 ) ), -1
                        x( ix ) = x( ix ) + temp*ap( k )
                        ix      = ix      - incx
   70                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*ap( kk - n + j )
                  end if
                  jx = jx - incx
                  kk = kk - ( n - j + 1 )
   80          continue
            end if
         end if
      else
c
c        form  x := a'*x  or  x := conjg( a' )*x.
c
         if( lsame( uplo, 'u' ) )then
            kk = ( n*( n + 1 ) )/2
            if( incx.eq.1 )then
               do 110, j = n, 1, -1
                  temp = x( j )
                  k    = kk     - 1
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*ap( kk )
                     do 90, i = j - 1, 1, -1
                        temp = temp + ap( k )*x( i )
                        k    = k    - 1
   90                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( ap( kk ) )
                     do 100, i = j - 1, 1, -1
                        temp = temp + conjg( ap( k ) )*x( i )
                        k    = k    - 1
  100                continue
                  end if
                  x( j ) = temp
                  kk     = kk   - j
  110          continue
            else
               jx = kx + ( n - 1 )*incx
               do 140, j = n, 1, -1
                  temp = x( jx )
                  ix   = jx
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*ap( kk )
                     do 120, k = kk - 1, kk - j + 1, -1
                        ix   = ix   - incx
                        temp = temp + ap( k )*x( ix )
  120                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( ap( kk ) )
                     do 130, k = kk - 1, kk - j + 1, -1
                        ix   = ix   - incx
                        temp = temp + conjg( ap( k ) )*x( ix )
  130                continue
                  end if
                  x( jx ) = temp
                  jx      = jx   - incx
                  kk      = kk   - j
  140          continue
            end if
         else
            kk = 1
            if( incx.eq.1 )then
               do 170, j = 1, n
                  temp = x( j )
                  k    = kk     + 1
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*ap( kk )
                     do 150, i = j + 1, n
                        temp = temp + ap( k )*x( i )
                        k    = k    + 1
  150                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( ap( kk ) )
                     do 160, i = j + 1, n
                        temp = temp + conjg( ap( k ) )*x( i )
                        k    = k    + 1
  160                continue
                  end if
                  x( j ) = temp
                  kk     = kk   + ( n - j + 1 )
  170          continue
            else
               jx = kx
               do 200, j = 1, n
                  temp = x( jx )
                  ix   = jx
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*ap( kk )
                     do 180, k = kk + 1, kk + n - j
                        ix   = ix   + incx
                        temp = temp + ap( k )*x( ix )
  180                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( ap( kk ) )
                     do 190, k = kk + 1, kk + n - j
                        ix   = ix   + incx
                        temp = temp + conjg( ap( k ) )*x( ix )
  190                continue
                  end if
                  x( jx ) = temp
                  jx      = jx   + incx
                  kk      = kk   + ( n - j + 1 )
  200          continue
            end if
         end if
      end if
c
      return
c
c     end of ctpmv .
c
      end
