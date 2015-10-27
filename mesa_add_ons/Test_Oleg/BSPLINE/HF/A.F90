!=======================================================================
  FUNCTION a(i,j,k)
!=======================================================================
!
!  determine the coefficient in the potential for electron i of
!  y^k(j,j)
!
!----------------------------------------------------------------------
!
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: i,j,k
     REAL(KIND=8) :: a

      if (i > nclosd .and. j > nclosd) then
         istart = ijptr(i-nclosd,j-nclosd) + 1
         a = coef(istart + k/2)
      else if (i == j) then
         c = sum(i) - 1.d0
         if (k == 0) then
            a = c
         else
            a = -c*ca(l(i),k)
         end if
      else if (k.eq.0) then
         a = sum(j)
      else
         a = 0.d0
      end if
  END FUNCTION a
