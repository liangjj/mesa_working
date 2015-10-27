      real*8           function slangt( norm, n, dl, d, du )
*
*  -- lapack auxiliary routine (version 1.1) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     february 29, 1992
*
*     .. scalar arguments ..
      character          norm
      integer            n
*     ..
*     .. array arguments ..
      real*8             d( * ), dl( * ), du( * )
*     ..
*
*  purpose
*  =======
*
*  slangt  returns the value of the one norm,  or the frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real tridiagonal matrix a.
*
*  description
*  ===========
*
*  slangt returns the value
*
*     slangt = ( max(abs(a(i,j))), norm = 'm' or 'm'
*              (
*              ( norm1(a),         norm = '1', 'o' or 'o'
*              (
*              ( normi(a),         norm = 'i' or 'i'
*              (
*              ( normf(a),         norm = 'f', 'f', 'e' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normi  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normf  denotes the  frobenius norm of a matrix (square root of sum of
*  squares).  note that  max(abs(a(i,j)))  is not a  matrix norm.
*
*  arguments
*  =========
*
*  norm    (input) character*1
*          specifies the value to be returned in slangt as described
*          above.
*
*  n       (input) integer
*          the order of the matrix a.  n >= 0.  when n = 0, slangt is
*          set to zero.
*
*  dl      (input) real array, dimension (n-1)
*          the (n-1) sub-diagonal elements of a.
*
*  d       (input) real array, dimension (n)
*          the diagonal elements of a.
*
*  du      (input) real array, dimension (n-1)
*          the (n-1) super-diagonal elements of a.
*
*  =====================================================================
*
*     .. parameters ..
      real*8             one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. local scalars ..
      integer            i
      real*8             anorm, scale, sum
*     ..
*     .. external functions ..
      logical            lsame
      external           lsame
*     ..
*     .. external subroutines ..
      external           slassq
*     ..
*     .. intrinsic functions ..
      intrinsic          abs, max, sqrt
*     ..
*     .. executable statements ..
*
      if( n.le.0 ) then
         anorm = zero
      else if( lsame( norm, 'm' ) ) then
*
*        find max(abs(a(i,j))).
*
         anorm = abs( d( n ) )
         do 10 i = 1, n - 1
            anorm = max( anorm, abs( dl( i ) ) )
            anorm = max( anorm, abs( d( i ) ) )
            anorm = max( anorm, abs( du( i ) ) )
   10    continue
      else if( lsame( norm, 'o' ) .or. norm.eq.'1' ) then
*
*        find norm1(a).
*
         if( n.eq.1 ) then
            anorm = abs( d( 1 ) )
         else
            anorm = max( abs( d( 1 ) )+abs( dl( 1 ) ),
     $              abs( d( n ) )+abs( du( n-1 ) ) )
            do 20 i = 2, n - 1
               anorm = max( anorm, abs( d( i ) )+abs( dl( i ) )+
     $                 abs( du( i-1 ) ) )
   20       continue
         end if
      else if( lsame( norm, 'i' ) ) then
*
*        find normi(a).
*
         if( n.eq.1 ) then
            anorm = abs( d( 1 ) )
         else
            anorm = max( abs( d( 1 ) )+abs( du( 1 ) ),
     $              abs( d( n ) )+abs( dl( n-1 ) ) )
            do 30 i = 2, n - 1
               anorm = max( anorm, abs( d( i ) )+abs( du( i ) )+
     $                 abs( dl( i-1 ) ) )
   30       continue
         end if
      else if( ( lsame( norm, 'f' ) ) .or. ( lsame( norm, 'e' ) ) ) then
*
*        find normf(a).
*
         scale = zero
         sum = one
         call slassq( n, d, 1, scale, sum )
         if( n.gt.1 ) then
            call slassq( n-1, dl, 1, scale, sum )
            call slassq( n-1, du, 1, scale, sum )
         end if
         anorm = scale*sqrt( sum )
      end if
*
      slangt = anorm
      return
*
*     end of slangt
*
      end
