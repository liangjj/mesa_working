!=======================================================================
   SUBROUTINE factrl(nfact)
!=======================================================================
!  Computes the log of the factorial function
!      gam(i) = log( gamma(i-1) ), where gamma(i) = factorial i-1
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nfact
    REAL(KIND=8), DIMENSION(*) :: gamma

      gamma = 1
      gam(1) = 0
      do  i=1,nfact-1
         gamma=i*gamma
         gam(i+1) = dlog(gamma)
      end do
      do i = nfact+1,(100)
         x = i-1
         gam(i) = gam(i-1) + dlog(x)
      end do
      
   END SUBROUTINE factrl
