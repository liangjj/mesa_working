      subroutine sgttrf( n, dl, d, du, du2, ipiv, det, info )
*
*  -- lapack routine (version 1.1) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     march 31, 1993 
*
*     .. scalar arguments ..
      integer            info, n
*     ..
*     .. array arguments ..
      integer            ipiv( * )
      real*8             d( * ), dl( * ), du( * ), du2( * ), det(2)
*     ..
*
*  purpose
*  =======
*
*  sgttrf computes an lu factorization of a real tridiagonal matrix a
*  using elimination with partial pivoting and row interchanges.
*  an option to include the calculation of the determinant is provided.
*
*  the factorization has the form
*     a = l * u
*  where l is a product of permutation and unit lower bidiagonal
*  matrices and u is upper triangular with nonzeros in only the main
*  diagonal and first two superdiagonals.
*
*  arguments
*  =========
*
*  n       (input) integer
*          the order of the matrix a.
*
*  dl      (input/output) real array, dimension (n-1)
*          on entry, dl must contain the (n-1) subdiagonal elements of
*          a.
*          on exit, dl is overwritten by the (n-1) multipliers that
*          define the matrix l from the lu factorization of a.
*
*  d       (input/output) real array, dimension (n)
*          on entry, d must contain the diagonal elements of a.
*          on exit, d is overwritten by the n diagonal elements of the
*          upper triangular matrix u from the lu factorization of a.
*
*  du      (input/output) real array, dimension (n-1)
*          on entry, du must contain the (n-1) superdiagonal elements
*          of a.
*          on exit, du is overwritten by the (n-1) elements of the first
*          superdiagonal of u.
*
*  du2     (output) real array, dimension (n-2)
*          on exit, du2 is overwritten by the (n-2) elements of the
*          second superdiagonal of u.
*
*  ipiv    (output) integer array, dimension (n)
*          the pivot indices; for 1 <= i <= n, row i of the matrix was
*          interchanged with row ipiv(i).  ipiv(i) will always be either
*          i or i+1; ipiv(i) = i indicates a row interchange was not
*          required.
*
*  info    (output) integer
*          = 0:  successful exit
*          < 0:  if info = -i, the i-th argument had an illegal value
*          > 0:  if info = i, u(i,i) is exactly zero. the factorization
*                has been completed, but the factor u is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. local scalars ..
      integer            i
      real*8             fact, temp
*
*     .. local logical ..
      logical            nodet
*     ..
*     .. intrinsic functions ..
      intrinsic          abs
*     ..
*     ..
*     .. parameters ..
      real*8             zero, one, ten, tol
      parameter          ( zero = 0.0d0, one=1.d0, ten=10.d0 )
      parameter          ( tol = 1.d-12 )
*     ..
*     .. executable statements ..
*
      if(info.eq.10) then
         nodet=.true.
      else
         nodet=.false.
      endif 
      info = 0
      if( n.lt.0 ) then
         info = -1
         call lnkerr( 'sgttrf called with n negative' )
      end if
*
*     quick return if possible
*
      if( n.eq.0 ) then
          call lnkerr( 'sgttrf called with n zero' )
      endif
*
*     initialize ipiv(i) = i
*
      do 10 i = 1, n
         ipiv( i ) = i
   10 continue
*
      do 20 i = 1, n - 1
         if( dl( i ).eq.zero ) then
*
*           subdiagonal is zero, no elimination is required.
*
            if( d( i ).eq.zero .and. info.eq.0 )
     $         info = i
            if( i.lt.n-1 )
     $         du2( i ) = zero
         else if( abs( d( i ) ).ge.abs( dl( i ) ) ) then
*
*           no row interchange required, eliminate dl(i)
*
            fact = dl( i ) / d( i )
            dl( i ) = fact
            d( i+1 ) = d( i+1 ) - fact*du( i )
            if( i.lt.n-1 )
     $         du2( i ) = zero
         else
*
*           interchange rows i and i+1, eliminate dl(i)
*
            fact = d( i ) / dl( i )
            d( i ) = dl( i )
            dl( i ) = fact
            temp = du( i )
            du( i ) = d( i+1 )
            d( i+1 ) = temp - fact*d( i+1 )
            if( i.lt.n-1 ) then
               du2( i ) = du( i+1 )
               du( i+1 ) = -fact*du( i+1 )
            end if
            ipiv( i ) = ipiv( i ) + 1
         end if
   20 continue
      if( d( n ).eq.zero .and. info.eq.0 ) then
         info = n
      endif
      if(nodet) then
         det(1) = one
         det(2) = zero
         do 50 i = 1, n                                                 
            if (ipiv(i) .ne. i) then
                det(1) = -det(1)                        
            endif
            det(1) = d(i)*det(1)
            if(abs(det(1)).le.tol) then
               det(1)=zero
            endif 
c           ...exit                                                        
            if (det(1) .ne. zero) then                             
                do while (abs(det(1)) .lt. one)                         
                    det(1) = ten*det(1)                                      
                    det(2) = det(2) - one
                enddo
                do while (abs(det(1)) .ge. ten)                           
                   det(1) = det(1)/ten                                      
                   det(2) = det(2) + one
                enddo
            else
                return                                                    
            endif
   50    continue

      endif
*
      return
*
*     end of sgttrf
*
      end
