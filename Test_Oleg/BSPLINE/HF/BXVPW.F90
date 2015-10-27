!=======================================================================
        subroutine bxvpw(c,b,v,y,w)
!=======================================================================
!   Computes y = c* b * v + y   where b is a symmetric, banded matrix
!   and v, w are vectors
!
!   Written by C. F. Fischer
!------------------------------------------------------------------------
!
!
!   on entry
!   --------
!     nt      the leading dimension of arrays.
!     k       the number of diagonals
!     n       the order of the matrix
!     c,ic    coefficients
!     b       the symmetric, banded matrix in column storage mode
!     v       vector
!     w       working array
!
!   on exit
!   -------
!     y       y = c*B*v +y
!-----------------------------------------------------------------------
!
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: i,j,k
     REAL(KIND=8), INTENT(IN), DIMENSION(nt,:) :: b
     REAL(KIND=8), INTENT(IN), DIMENSION(:) :: v,y
     REAL(KIND=8), INTENT(INOUT), DIMENSION(:) :: y

!
! ...   initialize the w array
!
	 w = 0.d0
!
!      .. contribution from sub-diagonals
!
       do i=1,k
         do j=k-i+1,n
	   w(j) = w(j)  +b(j,i)*v(j-k+i)
         end do
       end do
!
!      .. contribution from super-diagonals
!
       do i=1,k-1
         do j=1,n-k+i
	   w(j) = w(j) + b(j+k-i,i)*v(j+k-i)
         end do
       end do
!
       if ( c .ne. 1.d0) then
           y = y + c*w
       else
           y = y + w
       end if
    END SUBROUTINE bxvpw
