      subroutine sgttrs( trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb,
     $                   info )
*
*  -- lapack routine (version 1.1) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     march 31, 1993 
*
*     .. scalar arguments ..
      character          trans
      integer            info, ldb, n, nrhs
*     ..
*     .. array arguments ..
      integer            ipiv( * )
      real*8             b( ldb, * ), d( * ), dl( * ), du( * ), du2( * )
*     ..
*
*  purpose
*  =======
*
*  sgttrs solves one of the systems of equations
*     a*x = b  or  a'*x = b,
*  with a tridiagonal matrix a using the lu factorization computed
*  by sgttrf.
*
*  arguments
*  =========
*
*  trans   (input) character
*          specifies the form of the system of equations:
*          = 'n':  a * x = b  (no transpose)
*          = 't':  a'* x = b  (transpose)
*          = 'c':  a'* x = b  (conjugate transpose = transpose)
*
*  n       (input) integer
*          the order of the matrix a.
*
*  nrhs    (input) integer
*          the number of right hand sides, i.e., the number of columns
*          of the matrix b.  nrhs >= 0.
*
*  dl      (input) real array, dimension (n-1)
*          the (n-1) multipliers that define the matrix l from the
*          lu factorization of a.
*
*  d       (input) real array, dimension (n)
*          the n diagonal elements of the upper triangular matrix u from
*          the lu factorization of a.
*
*  du      (input) real array, dimension (n-1)
*          the (n-1) elements of the first superdiagonal of u.
*
*  du2     (input) real array, dimension (n-2)
*          the (n-2) elements of the second superdiagonal of u.
*
*  ipiv    (input) integer array, dimension (n)
*          the pivot indices; for 1 <= i <= n, row i of the matrix was
*          interchanged with row ipiv(i).  ipiv(i) will always be either
*          i or i+1; ipiv(i) = i indicates a row interchange was not
*          required.
*
*  b       (input/output) real array, dimension (ldb,nrhs)
*          on entry, the right hand side matrix b.
*          on exit, b is overwritten by the solution matrix x.
*
*  ldb     (input) integer
*          the leading dimension of the array b.  ldb >= max(1,n).
*
*  info    (output) integer
*          = 0:  successful exit
*          < 0:  if info = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. local scalars ..
      logical            notran
      integer            i, j
      real*8             temp
*     ..
*     .. external functions ..
      logical            lsame
      external           lsame
*     ..
*     .. external subroutines ..
      external           lnkerr
*     ..
*     .. intrinsic functions ..
      intrinsic          max
*     ..
*     .. executable statements ..
*
      info = 0
      notran = lsame( trans, 'n' )
      if( .not.notran .and. .not.lsame( trans, 't' ) .and. .not.
     $    lsame( trans, 'c' ) ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( nrhs.lt.0 ) then
         info = -3
      else if( ldb.lt.max( n, 1 ) ) then
         info = -10
      end if
      if( info.ne.0 ) then
         call lnkerr( 'sgttrs info .ne. zero' )
         return
      end if
*
*     quick return if possible
*
      if( n.eq.0 .or. nrhs.eq.0 )
     $   return
*
      if( notran ) then
*
*        solve a*x = b using the lu factorization of a,
*        overwriting each right hand side vector with its solution.
*
         do 30 j = 1, nrhs
*
*           solve l*x = b.
*
            do 10 i = 1, n - 1
               if( ipiv( i ).eq.i ) then
                  b( i+1, j ) = b( i+1, j ) - dl( i )*b( i, j )
               else
                  temp = b( i, j )
                  b( i, j ) = b( i+1, j )
                  b( i+1, j ) = temp - dl( i )*b( i, j )
               end if
   10       continue
*
*           solve u*x = b.
*
            b( n, j ) = b( n, j ) / d( n )
            if( n.gt.1 )
     $         b( n-1, j ) = ( b( n-1, j )-du( n-1 )*b( n, j ) ) /
     $                       d( n-1 )
            do 20 i = n - 2, 1, -1
               b( i, j ) = ( b( i, j )-du( i )*b( i+1, j )-du2( i )*
     $                     b( i+2, j ) ) / d( i )
   20       continue
   30    continue
      else
*
*        solve a' * x = b.
*
         do 60 j = 1, nrhs
*
*           solve u'*x = b.
*
            b( 1, j ) = b( 1, j ) / d( 1 )
            if( n.gt.1 )
     $         b( 2, j ) = ( b( 2, j )-du( 1 )*b( 1, j ) ) / d( 2 )
            do 40 i = 3, n
               b( i, j ) = ( b( i, j )-du( i-1 )*b( i-1, j )-du2( i-2 )*
     $                     b( i-2, j ) ) / d( i )
   40       continue
*
*           solve l'*x = b.
*
            do 50 i = n - 1, 1, -1
               if( ipiv( i ).eq.i ) then
                  b( i, j ) = b( i, j ) - dl( i )*b( i+1, j )
               else
                  temp = b( i+1, j )
                  b( i+1, j ) = b( i, j ) - dl( i )*temp
                  b( i, j ) = temp
               end if
   50       continue
   60    continue
      end if
*
*     end of sgttrs
*
      end
