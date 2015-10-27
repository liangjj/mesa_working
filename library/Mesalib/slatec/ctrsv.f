*deck ctrsv
      subroutine ctrsv (uplo, trans, diag, n, a, lda, x, incx)
c***begin prologue  ctrsv
c***purpose  solve a complex triangular system of equations.
c***library   slatec (blas)
c***category  d1b4
c***type      complex (strsv-s, dtrsv-d, ctrsv-c)
c***keywords  level 2 blas, linear algebra
c***author  dongarra, j. j., (anl)
c           du croz, j., (nag)
c           hammarling, s., (nag)
c           hanson, r. j., (snla)
c***description
c
c  ctrsv  solves one of the systems of equations
c
c     a*x = b,   or   a'*x = b,   or   conjg( a')*x = b,
c
c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular matrix.
c
c  no test for singularity or near-singularity is included in this
c  routine. such tests must be performed before calling this routine.
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
c           on entry, trans specifies the equations to be solved as
c           follows:
c
c              trans = 'n' or 'n'   a*x = b.
c
c              trans = 't' or 't'   a'*x = b.
c
c              trans = 'c' or 'c'   conjg( a' )*x = b.
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
c           element right-hand side vector b. on exit, x is overwritten
c           with the solution vector x.
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
c***end prologue  ctrsv
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
c***first executable statement  ctrsv
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
         call xerbla( 'ctrsv ', info )
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
c        form  x := inv( a )*x.
c
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/a( j, j )
                     temp = x( j )
                     do 10, i = j - 1, 1, -1
                        x( i ) = x( i ) - temp*a( i, j )
   10                continue
                  end if
   20          continue
            else
               jx = kx + ( n - 1 )*incx
               do 40, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/a( j, j )
                     temp = x( jx )
                     ix   = jx
                     do 30, i = j - 1, 1, -1
                        ix      = ix      - incx
                        x( ix ) = x( ix ) - temp*a( i, j )
   30                continue
                  end if
                  jx = jx - incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = 1, n
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/a( j, j )
                     temp = x( j )
                     do 50, i = j + 1, n
                        x( i ) = x( i ) - temp*a( i, j )
   50                continue
                  end if
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/a( j, j )
                     temp = x( jx )
                     ix   = jx
                     do 70, i = j + 1, n
                        ix      = ix      + incx
                        x( ix ) = x( ix ) - temp*a( i, j )
   70                continue
                  end if
                  jx = jx + incx
   80          continue
            end if
         end if
      else
c
c        form  x := inv( a' )*x  or  x := inv( conjg( a' ) )*x.
c
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 110, j = 1, n
                  temp = x( j )
                  if( noconj )then
                     do 90, i = 1, j - 1
                        temp = temp - a( i, j )*x( i )
   90                continue
                     if( nounit )
     $                  temp = temp/a( j, j )
                  else
                     do 100, i = 1, j - 1
                        temp = temp - conjg( a( i, j ) )*x( i )
  100                continue
                     if( nounit )
     $                  temp = temp/conjg( a( j, j ) )
                  end if
                  x( j ) = temp
  110          continue
            else
               jx = kx
               do 140, j = 1, n
                  ix   = kx
                  temp = x( jx )
                  if( noconj )then
                     do 120, i = 1, j - 1
                        temp = temp - a( i, j )*x( ix )
                        ix   = ix   + incx
  120                continue
                     if( nounit )
     $                  temp = temp/a( j, j )
                  else
                     do 130, i = 1, j - 1
                        temp = temp - conjg( a( i, j ) )*x( ix )
                        ix   = ix   + incx
  130                continue
                     if( nounit )
     $                  temp = temp/conjg( a( j, j ) )
                  end if
                  x( jx ) = temp
                  jx      = jx   + incx
  140          continue
            end if
         else
            if( incx.eq.1 )then
               do 170, j = n, 1, -1
                  temp = x( j )
                  if( noconj )then
                     do 150, i = n, j + 1, -1
                        temp = temp - a( i, j )*x( i )
  150                continue
                     if( nounit )
     $                  temp = temp/a( j, j )
                  else
                     do 160, i = n, j + 1, -1
                        temp = temp - conjg( a( i, j ) )*x( i )
  160                continue
                     if( nounit )
     $                  temp = temp/conjg( a( j, j ) )
                  end if
                  x( j ) = temp
  170          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 200, j = n, 1, -1
                  ix   = kx
                  temp = x( jx )
                  if( noconj )then
                     do 180, i = n, j + 1, -1
                        temp = temp - a( i, j )*x( ix )
                        ix   = ix   - incx
  180                continue
                     if( nounit )
     $                  temp = temp/a( j, j )
                  else
                     do 190, i = n, j + 1, -1
                        temp = temp - conjg( a( i, j ) )*x( ix )
                        ix   = ix   - incx
  190                continue
                     if( nounit )
     $                  temp = temp/conjg( a( j, j ) )
                  end if
                  x( jx ) = temp
                  jx      = jx   - incx
  200          continue
            end if
         end if
      end if
c
      return
c
c     end of ctrsv .
c
      end
