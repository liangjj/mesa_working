*deck cgerc
      subroutine cgerc (m, n, alpha, x, incx, y, incy, a, lda)
c***begin prologue  cgerc
c***purpose  perform conjugated rank 1 update of a complex general
c            matrix.
c***library   slatec (blas)
c***category  d1b4
c***type      complex (sgerc-s, dgerc-d, cgerc-c)
c***keywords  level 2 blas, linear algebra
c***author  dongarra, j. j., (anl)
c           du croz, j., (nag)
c           hammarling, s., (nag)
c           hanson, r. j., (snla)
c***description
c
c  cgerc  performs the rank 1 operation
c
c     a := alpha*x*conjg( y') + a,
c
c  where alpha is a scalar, x is an m element vector, y is an n element
c  vector and a is an m by n matrix.
c
c  parameters
c  ==========
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
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( m - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the m
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  y      - complex          array of dimension at least
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
c  a      - complex          array of dimension ( lda, n ).
c           before entry, the leading m by n part of the array a must
c           contain the matrix of coefficients. on exit, a is
c           overwritten by the updated matrix.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, m ).
c           unchanged on exit.
c
c***references  dongarra, j. j., du croz, j., hammarling, s., and
c                 hanson, r. j.  an extended set of fortran basic linear
c                 algebra subprograms.  acm toms, vol. 14, no. 1,
c                 pp. 1-17, march 1988.
c***routines called  xerbla
c***revision history  (yymmdd)
c   861022  date written
c   910605  modified to meet slatec prologue standards.  only comment
c           lines were modified.  (bks)
c***end prologue  cgerc
c     .. scalar arguments ..
      complex            alpha
      integer            incx, incy, lda, m, n
c     .. array arguments ..
      complex            a( lda, * ), x( * ), y( * )
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jy, kx
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max
c***first executable statement  cgerc
c
c     test the input parameters.
c
      info = 0
      if     ( m.lt.0 )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( incy.eq.0 )then
         info = 7
      else if( lda.lt.max( 1, m ) )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'cgerc ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
      if( incy.gt.0 )then
         jy = 1
      else
         jy = 1 - ( n - 1 )*incy
      end if
      if( incx.eq.1 )then
         do 20, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*conjg( y( jy ) )
               do 10, i = 1, m
                  a( i, j ) = a( i, j ) + x( i )*temp
   10          continue
            end if
            jy = jy + incy
   20    continue
      else
         if( incx.gt.0 )then
            kx = 1
         else
            kx = 1 - ( m - 1 )*incx
         end if
         do 40, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*conjg( y( jy ) )
               ix   = kx
               do 30, i = 1, m
                  a( i, j ) = a( i, j ) + x( ix )*temp
                  ix        = ix        + incx
   30          continue
            end if
            jy = jy + incy
   40    continue
      end if
c
      return
c
c     end of cgerc .
c
      end
