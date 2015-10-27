!deck Re_Poly.f
!***begin prologue     Re_Poly
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
!***end prologue       Re_Poly

  SUBROUTINE Re_Poly(p,dp,ddp,inv_sqrt_wt,n)
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)          :: p
  REAL*8, DIMENSION(:,:)          :: dp
  REAL*8, DIMENSION(:,:)          :: ddp
  REAL*8, DIMENSION (:)           :: inv_sqrt_wt
  INTEGER                         :: n
  INTEGER                         :: i
!
!
  DO i = 1, n
       p(:, i ) =    p(:, i )  * inv_sqrt_wt(i)
      dp(:, i ) =   dp(:, i )  * inv_sqrt_wt(i)
     ddp(:, i ) =  ddp(:, i )  * inv_sqrt_wt(i)
  END DO
END SUBROUTINE Re_Poly
