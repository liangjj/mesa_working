
!=======================================================================
   FUNCTION cb(l,lp,k) 
!=======================================================================
!  Compute the average exchange interaction
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN): l,lp,k
    REAL(KIND=8):: cb

      cb = rme(l,lp,k)**2/(2*(2*l+1)*(2*lp+1))

    END FUNCTION cb
