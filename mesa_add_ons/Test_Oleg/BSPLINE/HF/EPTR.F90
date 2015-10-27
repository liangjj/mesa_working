!=======================================================================
  SUBROUTINE eptr(el, elsymb, iel) 
!=======================================================================
!
!   Determines the position of the electron in the electron list
!   Zero if not found.
!      
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    CHARACTER(LEN=*), DIMENSION(*), INTENT(IN) :: el
    CHARACTER(LEN=*), INTENT(IN) :: elsymb
    INTEGER, INTENT(OUT) :: iel

    INTEGER :: n, i
    n = size(el)
    iel = 0
    do  i=1,n
      if (el(i) .eq. elsymb ) then
        iel = i
        return
      endif
    end do
    end
