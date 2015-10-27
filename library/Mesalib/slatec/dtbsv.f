*deck dtbsv
      subroutine dtbsv (uplo, trans, diag, n, k, a, lda, x, incx)
c***begin prologue  dtbsv
c***purpose  solve one of the systems of equations.
c***library   slatec (blas)
c***category  d1b4
c***type      double precision (stbsv-s, dtbsv-d, ctbsv-c)
c***keywords  level 2 blas, linear algebra
c***author  dongarra, j. j., (anl)
c           du croz, j., (nag)
c           hammarling, s., (nag)
c           hanson, r. j., (snla)
c***description
c
c  dtbsv  solves one of the systems of equations
c
c     a*x = b,   or   a'*x = b,
c
c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular band matrix, with ( k + 1)
c  diagonals.
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
c              trans = 'c' or 'c'   a'*x = b.
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
c  k      - integer.
c           on entry with uplo = 'u' or 'u', k specifies the number of
c           super-diagonals of the matrix a.
c           on entry with uplo = 'l' or 'l', k specifies the number of
c           sub-diagonals of the matrix a.
c           k must satisfy  0 .le. k.
c           unchanged on exit.
c
c  a      - double precision array of dimension ( lda, n ).
c           before entry with uplo = 'u' or 'u', the leading ( k + 1 )
c           by n part of the array a must contain the upper triangular
c           band part of the matrix of coefficients, supplied column by
c           column, with the leading diagonal of the matrix in row
c           ( k + 1 ) of the array, the first super-diagonal starting at
c           position 2 in row k, and so on. the top left k by k triangle
c           of the array a is not referenced.
c           the following program segment will transfer an upper
c           triangular band matrix from conventional full matrix storage
c           to band storage:
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
c           band part of the matrix of coefficients, supplied column by
c           column, with the leading diagonal of the matrix in row 1 of
c           the array, the first sub-diagonal starting at position 1 in
c           row 2, and so on. the bottom right k by k triangle of the
c           array a is not referenced.
c           the following program segment will transfer a lower
c           triangular band matrix from conventional full matrix storage
c           to band storage:
c
c                 do 20, j = 1, n
c                    m = 1 - j
c                    do 10, i = j, min( n, j + k )
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           note that when diag = 'u' or 'u' the elements of the array a
c           corresponding to the diagonal elements of the matrix are not
c           referenced, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           ( k + 1 ).
c           unchanged on exit.
c
c  x      - double precision array of dimension at least
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
c***end prologue  dtbsv
c     .. scalar arguments ..
      integer            incx, k, lda, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      double precision   a( lda, * ), x( * )
c     .. parameters ..
      double precision   zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      double precision   temp
      integer            i, info, ix, j, jx, kplus1, kx, l
      logical            nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max, min
c***first executable statement  dtbsv
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
      else if( k.lt.0 )then
         info = 5
      else if( lda.lt.( k + 1 ) )then
         info = 7
      else if( incx.eq.0 )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'dtbsv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      nounit = lsame( diag, 'n' )
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
c     accessed by sequentially with one pass through a.
c
      if( lsame( trans, 'n' ) )then
c
c        form  x := inv( a )*x.
c
         if( lsame( uplo, 'u' ) )then
            kplus1 = k + 1
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     l = kplus1 - j
                     if( nounit )
     $                  x( j ) = x( j )/a( kplus1, j )
                     temp = x( j )
                     do 10, i = j - 1, max( 1, j - k ), -1
                        x( i ) = x( i ) - temp*a( l + i, j )
   10                continue
                  end if
   20          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 40, j = n, 1, -1
                  kx = kx - incx
                  if( x( jx ).ne.zero )then
                     ix = kx
                     l  = kplus1 - j
                     if( nounit )
     $                  x( jx ) = x( jx )/a( kplus1, j )
                     temp = x( jx )
                     do 30, i = j - 1, max( 1, j - k ), -1
                        x( ix ) = x( ix ) - temp*a( l + i, j )
                        ix      = ix      - incx
   30                continue
                  end if
                  jx = jx - incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = 1, n
                  if( x( j ).ne.zero )then
                     l = 1 - j
                     if( nounit )
     $                  x( j ) = x( j )/a( 1, j )
                     temp = x( j )
                     do 50, i = j + 1, min( n, j + k )
                        x( i ) = x( i ) - temp*a( l + i, j )
   50                continue
                  end if
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  kx = kx + incx
                  if( x( jx ).ne.zero )then
                     ix = kx
                     l  = 1  - j
                     if( nounit )
     $                  x( jx ) = x( jx )/a( 1, j )
                     temp = x( jx )
                     do 70, i = j + 1, min( n, j + k )
                        x( ix ) = x( ix ) - temp*a( l + i, j )
                        ix      = ix      + incx
   70                continue
                  end if
                  jx = jx + incx
   80          continue
            end if
         end if
      else
c
c        form  x := inv( a')*x.
c
         if( lsame( uplo, 'u' ) )then
            kplus1 = k + 1
            if( incx.eq.1 )then
               do 100, j = 1, n
                  temp = x( j )
                  l    = kplus1 - j
                  do 90, i = max( 1, j - k ), j - 1
                     temp = temp - a( l + i, j )*x( i )
   90             continue
                  if( nounit )
     $               temp = temp/a( kplus1, j )
                  x( j ) = temp
  100          continue
            else
               jx = kx
               do 120, j = 1, n
                  temp = x( jx )
                  ix   = kx
                  l    = kplus1  - j
                  do 110, i = max( 1, j - k ), j - 1
                     temp = temp - a( l + i, j )*x( ix )
                     ix   = ix   + incx
  110             continue
                  if( nounit )
     $               temp = temp/a( kplus1, j )
                  x( jx ) = temp
                  jx      = jx   + incx
                  if( j.gt.k )
     $               kx = kx + incx
  120          continue
            end if
         else
            if( incx.eq.1 )then
               do 140, j = n, 1, -1
                  temp = x( j )
                  l    = 1      - j
                  do 130, i = min( n, j + k ), j + 1, -1
                     temp = temp - a( l + i, j )*x( i )
  130             continue
                  if( nounit )
     $               temp = temp/a( 1, j )
                  x( j ) = temp
  140          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 160, j = n, 1, -1
                  temp = x( jx )
                  ix   = kx
                  l    = 1       - j
                  do 150, i = min( n, j + k ), j + 1, -1
                     temp = temp - a( l + i, j )*x( ix )
                     ix   = ix   - incx
  150             continue
                  if( nounit )
     $               temp = temp/a( 1, j )
                  x( jx ) = temp
                  jx      = jx   - incx
                  if( ( n - j ).ge.k )
     $               kx = kx - incx
  160          continue
            end if
         end if
      end if
c
      return
c
c     end of dtbsv .
c
      end
