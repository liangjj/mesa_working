      integer function idmax(n,dx,incx)
!=======================================================================
  FUNCTION idmax(n,dx,incx) 
!=======================================================================
!     finds the index of element having largest value.
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) ::n, incx
    REAL(KIND=8), INTENT(IN) :: dx

    INTEGER :: idmax, i

    if( n < 1 ) then
      idmax = 0
    else if (n == 1) then
      idmax = 1
    else
      if (incx > 1) then
        ! code for increment not equal to 1
        ix = 1
        dmax = dx(1)
        ix = ix + incx
        do i = 2,n
         if(dx(ix) > dmax) then
           idmax = i
           dmax = dx(ix)
          end if
       	end do
      else
        !code for increment equal to 1
	idmax = 1
        dmax = dx(1)
        do i = 2,n
         if(dx(i) > dmax) then
           idmax = i
           dmax = dx(i)
	 end if
	end do
      end if
   END FUNCTION idmax
