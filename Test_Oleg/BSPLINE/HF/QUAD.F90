!=======================================================================
  FUNCTION quad(a,b) 
!=======================================================================
!       Returns the value of <a, S b> where a and b are spline
!     expansion coefficients
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(ns):: a, b

      quad = 0
      do i = 2,ns
	quad = quad + a(i)*b(i)*sb(i,ks)
      end do
      do m = 2,ks
        ! i-m+1  < m
	do i = m+1,ns
	  quad = quad +(a(i)*b(i-m+1)+a(i-m+1)*b(i))*sb(i,ks-m+1)
        end do
      end do
    END FUNCTION quad
