      subroutine cgetrs( trans, n, nrhs, a, lda, ipiv, b, ldb, info )
*
*  -- lapack routine (version 1.1) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     march 31, 1993
*
*     .. scalar arguments ..
      character          trans
      integer            info, lda, ldb, n, nrhs
*     ..
*     .. array arguments ..
      integer            ipiv( * )
      complex*16         a( lda, * ), b( ldb, * )
*     ..
*
*  purpose
*  =======
*
*  cgetrs solves a system of linear equations
*     a * x = b,  a**t * x = b,  or  a**h * x = b
*  with a general n-by-n matrix a using the lu factorization computed
*  by cgetrf.
*
*  arguments
*  =========
*
*  trans   (input) character*1
*          specifies the form of the system of equations:
*          = 'n':  a * x = b     (no transpose)
*          = 't':  a**t * x = b  (transpose)
*          = 'c':  a**h * x = b  (conjugate transpose)
*
*  n       (input) integer
*          the order of the matrix a.  n >= 0.
*
*  nrhs    (input) integer
*          the number of right hand sides, i.e., the number of columns
*          of the matrix b.  nrhs >= 0.
*
*  a       (input) complex*16 array, dimension (lda,n)
*          the factors l and u from the factorization a = p*l*u
*          as computed by zgetrf.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(1,n).
*
*  ipiv    (input) integer array, dimension (n)
*          the pivot indices from zgetrf; for 1<=i<=n, row i of the
*          matrix was interchanged with row ipiv(i).
*
*  b       (input/output) complex*16 array, dimension (ldb,nrhs)
*          on entry, the right hand side matrix b.
*          on exit, the solution matrix x.
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
*     .. parameters ..
      complex*16         one
      parameter          ( one = 1.0d+0 )
*     ..
*     .. local scalars ..
      logical            notran
*     ..
*     .. external functions ..
      logical            lsame
      external           lsame
*     ..
*     .. external subroutines ..
      external           xerbla, zlaswp, ztrsm
*     ..
*     .. intrinsic functions ..
      intrinsic          max
*     ..
*     .. executable statements ..
*
#ifdef LOCALLAPACK
*     test the input parameters.
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
      else if( lda.lt.max( 1, n ) ) then
         info = -5
      else if( ldb.lt.max( 1, n ) ) then
         info = -8
      end if
      if( info.ne.0 ) then
         call xerbla( 'zgetrs', -info )
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
*        solve a * x = b.
*
*        apply row interchanges to the right hand sides.
*
         call zlaswp( nrhs, b, ldb, 1, n, ipiv, 1 )
*
*        solve l*x = b, overwriting b with x.
*
         call ztrsm( 'left', 'lower', 'no transpose', 'unit', n, nrhs,
     $               one, a, lda, b, ldb )
*
*        solve u*x = b, overwriting b with x.
*
         call ztrsm( 'left', 'upper', 'no transpose', 'non-unit', n,
     $               nrhs, one, a, lda, b, ldb )
      else
*
*        solve a**t * x = b  or a**h * x = b.
*
*        solve u'*x = b, overwriting b with x.
*
         call ztrsm( 'left', 'upper', trans, 'non-unit', n, nrhs, one,
     $               a, lda, b, ldb )
*
*        solve l'*x = b, overwriting b with x.
*
         call ztrsm( 'left', 'lower', trans, 'unit', n, nrhs, one, a,
     $               lda, b, ldb )
*
*        apply row interchanges to the solution vectors.
*
         call zlaswp( nrhs, b, ldb, 1, n, ipiv, -1 )
      end if
*
      return
*
*     end of zgetrs
*
      end
#endif LOCALLAPACK
#ifdef MACHINELAPACK
      call zgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info )
      return
      end
#endif MACHINELAPACK
