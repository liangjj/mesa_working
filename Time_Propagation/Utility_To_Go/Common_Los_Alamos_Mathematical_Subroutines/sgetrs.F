      subroutine sgetrs( trans, n, nrhs, a, lda, ipiv, b, ldb, info )
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
      double precision   a( lda, * ), b( ldb, * )
*     ..
*
*  purpose
*  =======
*
*  dgetrs solves a system of linear equations
*     a * x = b  or  a' * x = b
*  with a general n-by-n matrix a using the lu factorization computed
*  by sgetrf.
*
*  arguments
*  =========
*
*  trans   (input) character*1
*          specifies the form of the system of equations:
*          = 'n':  a * x = b  (no transpose)
*          = 't':  a'* x = b  (transpose)
*          = 'c':  a'* x = b  (conjugate transpose = transpose)
*
*  n       (input) integer
*          the order of the matrix a.  n >= 0.
*
*  nrhs    (input) integer
*          the number of right hand sides, i.e., the number of columns
*          of the matrix b.  nrhs >= 0.
*
*  a       (input) double precision array, dimension (lda,n)
*          the factors l and u from the factorization a = p*l*u
*          as computed by dgetrf.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(1,n).
*
*  ipiv    (input) integer array, dimension (n)
*          the pivot indices from dgetrf; for 1<=i<=n, row i of the
*          matrix was interchanged with row ipiv(i).
*
*  b       (input/output) double precision array, dimension (ldb,nrhs)
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
      double precision   one
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
      external           dlaswp, dtrsm, xerbla
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
         call xerbla( 'dgetrs', -info )
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
         call dlaswp( nrhs, b, ldb, 1, n, ipiv, 1 )
*
*        solve l*x = b, overwriting b with x.
*
         call dtrsm( 'left', 'lower', 'no transpose', 'unit', n, nrhs,
     $               one, a, lda, b, ldb )
*
*        solve u*x = b, overwriting b with x.
*
         call dtrsm( 'left', 'upper', 'no transpose', 'non-unit', n,
     $               nrhs, one, a, lda, b, ldb )
      else
*
*        solve a' * x = b.
*
*        solve u'*x = b, overwriting b with x.
*
         call dtrsm( 'left', 'upper', 'transpose', 'non-unit', n, nrhs,
     $               one, a, lda, b, ldb )
*
*        solve l'*x = b, overwriting b with x.
*
         call dtrsm( 'left', 'lower', 'transpose', 'unit', n, nrhs, one,
     $               a, lda, b, ldb )
*
*        apply row interchanges to the solution vectors.
*
         call dlaswp( nrhs, b, ldb, 1, n, ipiv, -1 )
      end if
      return
      end
#endif LOCALLAPACK
#ifdef MACHINELAPACK
      call dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info )
      return
      end
#endif MACHINELAPACK
