!=======================================================================
  FUNCTION b(i,j,k) 
!=======================================================================
!
!   determine the coefficient of the y^k(i,j)p(j) term in the exchange
!   expression of electron i
!
!----------------------------------------------------------------------
!
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: i,j,k
   REAL(KIND=8) :: b

 
   if (i == j) then
       b = 0.d0
   else if (i > nclosd .and. j > nclosd) then
 
    !  ll is the number of direct terms
    !  istart the beginning of the exchange terms
 
    ll = min(l(i),l(j)) + 1
    istart = ijptr(i-nclosd,j-nclosd) + 1 + ll
    kk = (k - abs(l(i)-l(j)))/2
    b = coef(istart + kk)
  else
    b = -sum(j)*cb(l(i),l(j),k)
  end if
  END FUNCTION b
