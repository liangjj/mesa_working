!deck ReGrid.f
!***begin prologue     ReGrid
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Modify the weights at the ends of the interval to
!***                   reflect that there are bridge functions present. 
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       ReGrid

  SUBROUTINE ReGrid(q,wt,inv_sqrt_wt,wt_left_end,wt_right_end,n,region)
  USE dvr_global,            ONLY : nreg
  IMPLICIT NONE
  REAL*8, DIMENSION (:)             :: q
  REAL*8, DIMENSION (:)             :: wt
  REAL*8, DIMENSION (:)             :: inv_sqrt_wt
  REAL*8                            :: wt_end
  INTEGER                           :: n
  INTEGER                           :: region
!
!
  IF ( nreg == 1) THEN
       i = 1
       inv_sqrt_wt(1 : n) = Sqrt ( 1.d0 / reg_grid(i)%wt( 1 : n) )
       Return
  END IF
  IF ( region == 1) THEN
!
!      Modify the last weight and then get the inverse square roots.
!
       wt( n ) = wt( n ) + wt_right_end
       inv_sqrt_wt( 1 : n ) = 1.d0 / sqrt ( wt( 1 : n) )     
  ELSE IF ( region == nreg ) THEN
!
!      
!
       wt( 1 ) = wt_left_end
       inv_sqrt_wt( 1 : n ) = 1.d0 / sqrt ( wt( 1 : n ) )
  ELSE
!      
!      Modify the last weight and then get the inverse square roots.
!
       wt( 1 )  = wt_left_end
       wt( n )  = wt( n ) + wt_right_end
       inv_sqrt_wt( 1 : n )  =  1.d0 / sqrt ( wt( 1 : n ) )
  END IF
END SUBROUTINE ReGrid
