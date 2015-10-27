!deck Renorm.f
!***begin prologue     Renorm
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            1. calculate piecewise lobatto dvr functions and
!***                      their one-body matrices
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Renorm

  SUBROUTINE Renorm_Eta( grid, region, n)
  IMPLICIT NONE
  TYPE(coordinates)         :: grid
  INTEGER                   :: region
  REAL*8                    :: one = 1.d0
  REAL*8                    :: fac
  INTEGER                   :: i
!
!
  DO i = 1, n
     fac = Sqrt ( one / ( one - grid%pt_wt(region)%q(i) * grid%pt_wt(region)%q(i) * q(i) ) )
         p_o(:,i)  = fac * p_e(:,j) 
         dp_o(:,i) = fac * dp_e(:,i) 
         ddp_o(:,i) = fac * ddp_e(:,i) 
      END DO
  ELSE IF (coord == 'xi' ) THEN
      DO i= 1, n
         fac = Sqrt ( one / ( q(i) * q(i) - one ) )
         p_o(:,i) = fac * p_e(:,i) 
         dp_o(:,i) = fac * dp_e(:,i) 
         ddpr_o(:,i) = fac * ddp_e(:,i) 
      END DO
  END IF
END SUBROUTINE Renorm_Eta
