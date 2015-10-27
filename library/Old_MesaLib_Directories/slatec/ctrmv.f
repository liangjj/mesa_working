*deck ctrmv
      subroutine ctrmv (uplo, trans, diag, n, a, lda, x, incx)
c***begin prologue  ctrmv
c***purpose  multiply a complex vector by a complex triangular matrix.
c***library   slatec (blas)
c***category  d1b4
c***type      complex (strmv-s, dtrmv-d, ctrmv-c)
c***keywords  level 2 blas, linear algebra
c***author  dongarra, j. j., (anl)
c           du croz, j., (nag)
c           hammarling, s., (nag)
c           hanson, r. j., (snla)
c***description
c
c  ctrmv  performs one of the matrix-vector operations
c
c     x := a*x,   or   x := a'*x,   or   x := conjg( a')*x,
c
c  where x is an n element vector and  a is an n by n unit, or non-unit,
c  upper or lower triangular matrix.
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
c  a      - complex          array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular matrix and the strictly lower triangular part of
c           a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular matrix and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced either, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
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
c***end prologue  ctrmv
c     .. scalar arguments ..
      integer            incx, lda, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      complex            a( lda, * ), x( * )
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jx, kx
      logical            noconj, nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max
c***first executable statement  ctrmv
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
      else if( lda.lt.max( 1, n ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      end if
      if( info.ne.0 )then
         call xerbla( 'ctrmv ', info )
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
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
      if( lsame( trans, 'n' ) )then
c
c        form  x := a*x.
c
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 20, j = 1, n
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     do 10, i = 1, j - 1
                        x( i ) = x( i ) + temp*a( i, j )
   10                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( j, j )
                  end if
   20          continue
            else
               jx = kx
               do 40, j = 1, n
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 30, i = 1, j - 1
                        x( ix ) = x( ix ) + temp*a( i, j )
                        ix      = ix      + incx
   30                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( j, j )
                  end if
                  jx = jx + incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     do 50, i = n, j + 1, -1
                        x( i ) = x( i ) + temp*a( i, j )
   50                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( j, j )
                  end if
   60          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 80, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 70, i = n, j + 1, -1
                        x( ix ) = x( ix ) + temp*a( i, j )
                        ix      = ix      - incx
   70                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( j, j )
                  end if
                  jx = jx - incx
   80          continue
            end if
         end if
      else
c
c        form  x := a'*x  or  x := conjg( a' )*x.
c
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 110, j = n, 1, -1
                  temp = x( j )
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*a( j, j )
                     do 90, i = j - 1, 1, -1
                        temp = temp + a( i, j )*x( i )
   90                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( a( j, j ) )
                     do 100, i = j - 1, 1, -1
                        temp = temp + conjg( a( i, j ) )*x( i )
  100                continue
                  end if
                  x( j ) = temp
  110          continue
            else
               jx = kx + ( n - 1 )*incx
               do 140, j = n, 1, -1
                  temp = x( jx )
                  ix   = jx
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*a( j, j )
                     do 120, i = j - 1, 1, -1
                        ix   = ix   - incx
                        temp = temp + a( i, j )*x( ix )
  120                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( a( j, j ) )
                     do 130, i = j - 1, 1, -1
                        ix   = ix   - incx
                        temp = temp + conjg( a( i, j ) )*x( ix )
  130                continue
                  end if
                  x( jx ) = temp
                  jx      = jx   - incx
  140          continue
            end if
         else
            if( incx.eq.1 )then
               do 170, j = 1, n
                  temp = x( j )
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*a( j, j )
                     do 150, i = j + 1, n
                        temp = temp + a( i, j )*x( i )
  150                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( a( j, j ) )
                     do 160, i = j + 1, n
                        temp = temp + conjg( a( i, j ) )*x( i )
  160                continue
                  end if
                  x( j ) = temp
  170          continue
            else
               jx = kx
               do 200, j = 1, n
                  temp = x( jx )
                  ix   = jx
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*a( j, j )
                     do 180, i = j + 1, n
                        ix   = ix   + incx
                        temp = temp + a( i, j )*x( ix )
  180                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( a( j, j ) )
                     do 190, i = j + 1, n
                        ix   = ix   + incx
                        temp = temp + conjg( a( i, j ) )*x( ix )
  190                continue
                  end if
                  x( jx ) = temp
                  jx      = jx   + incx
  200          continue
            end if
         end if
      end if
c
      return
c
c     end of ctrmv .
c
      end
