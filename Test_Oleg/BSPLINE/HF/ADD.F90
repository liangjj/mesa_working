!=======================================================================
  SUBROUTINE add(c,k,i,j,first) 
!=======================================================================
!   add a slater integral to the data structure associated with the
!   energy expression
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k,i,j
    LOGICAL, INTENT(IN) :: first
    REAL(KIND=8), INTENT(IN) :: c

    ip = ijptr(i-nclosd,j-nclosd)
   
    if (first) then
       coef(ip+k/2+1) = c/sum(i) + coef(ip+k/2+1)
    else
       ip = ip + min(l(i),l(j)) +1 + (k-abs(l(i)-l(j)))/2 + 1
       coef(ip) = coef(ip) + c/sum(i)
    end if

  END FUNCTION add
