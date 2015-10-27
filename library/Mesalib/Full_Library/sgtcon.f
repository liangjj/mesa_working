      subroutine sgtcon( norm, n, dl, d, du, du2, ipiv, anorm, rcond,
     $                   work, iwork, info )
*
*  -- lapack routine (version 1.1) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     march 31, 1993 
*
*     .. scalar arguments ..
      character          norm
      integer            info, n
      real*8             anorm, rcond
*     ..
*     .. array arguments ..
      integer            ipiv( * ), iwork( * )
      real*8             d( * ), dl( * ), du( * ), du2( * ), work( * )
*     ..
*
*  purpose
*  =======
*
*  sgtcon estimates the reciprocal of the condition number of a real
*  tridiagonal matrix a using the lu factorization as computed by
*  sgttrf.
*
*  an estimate is obtained for norm(inv(a)), and the reciprocal of the
*  condition number is computed as rcond = 1 / (anorm * norm(inv(a))).
*
*  arguments
*  =========
*
*  norm    (input) character*1
*          specifies whether the 1-norm condition number or the
*          infinity-norm condition number is required:
*          = '1' or 'o':  1-norm;
*          = 'i':         infinity-norm.
*
*  n       (input) integer
*          the order of the matrix a.  n >= 0.
*
*  dl      (input) real array, dimension (n-1)
*          the (n-1) multipliers that define the matrix l from the
*          lu factorization of a as computed by sgttrf.
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
*  anorm   (input) real
*          the 1-norm of the original matrix a.
*
*  rcond   (output) real
*          the reciprocal of the condition number of the matrix a,
*          computed as rcond = 1/(anorm * ainvnm), where ainvnm is an
*          estimate of the 1-norm of inv(a) computed in this routine.
*
*  work    (workspace) real array, dimension (2*n)
*
*  iwork   (workspace) integer array, dimension (n)
*
*  info    (output) integer
*          = 0:  successful exit
*          < 0:  if info = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. parameters ..
      real*8             one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. local scalars ..
      logical            onenrm
      integer            i, kase, kase1
      real*8             ainvnm
*     ..
*     .. external functions ..
      logical            lsame
      external           lsame
*     ..
*     .. external subroutines ..
      external           sgttrs, slacon, lnkerr
*     ..
*     .. executable statements ..
*
*     test the input arguments.
*
      info = 0
      onenrm = norm.eq.'1' .or. lsame( norm, 'o' )
      if( .not.onenrm .and. .not.lsame( norm, 'i' ) ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( anorm.lt.zero ) then
         info = -8
      end if
      if( info.ne.0 ) then
         call lnkerr( 'sgtcon info .ne. zero' )
         return
      end if
*
*     quick return if possible
*
      rcond = zero
      if( n.eq.0 ) then
         rcond = one
         return
      else if( anorm.eq.zero ) then
         return
      end if
*
*     check that d(1:n) is non-zero.
*
      do 10 i = 1, n
         if( d( i ).eq.zero )
     $      return
   10 continue
*
      ainvnm = zero
      if( onenrm ) then
         kase1 = 1
      else
         kase1 = 2
      end if
      kase = 0
   20 continue
      call slacon( n, work( n+1 ), work, iwork, ainvnm, kase )
      if( kase.ne.0 ) then
         if( kase.eq.kase1 ) then
*
*           multiply by inv(u)*inv(l).
*
            call sgttrs( 'no transpose', n, 1, dl, d, du, du2, ipiv,
     $                   work, n, info )
         else
*
*           multiply by inv(l')*inv(u').
*
            call sgttrs( 'transpose', n, 1, dl, d, du, du2, ipiv, work,
     $                   n, info )
         end if
         go to 20
      end if
*
*     compute the estimate of the reciprocal condition number.
*
      if( ainvnm.ne.zero )
     $   rcond = ( one / ainvnm ) / anorm
*
      return
*
*     end of sgtcon
*
      end
