      
!=======================================================================
   FUNCTION ca(l,k) 
!=======================================================================
! Computes the direct average interaction
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: l, k
    REAL(KIND=8) :: ca

      ca = rme(l,l,k)**2

    END FUNCTION ca
