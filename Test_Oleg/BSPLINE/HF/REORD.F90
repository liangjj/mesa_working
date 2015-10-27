
!=======================================================================
  SUBROUTINE reord(el, elc, nwf, ierr) 
!=======================================================================
!
!   Reorder the list of first appearance so that the FUNCTIONs (orbitals)
!   to be varied appear last in the list.
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    CHARACTER(LEN=*), DIMENSION(nwf), INTENT(INOUT) :: el
    CHARACTER(LEN=*), INTENT(IN) :: elc 
*
    call eptr(el, elc, i)
    if (i <> 0) then
      do j = i, nwf-1
        el(j) = el(j+1)
      end do
      el(nwf) = elc
      ierr = 0
    else
      ierr = 1
    end if
  END SUBROUTINE reord


