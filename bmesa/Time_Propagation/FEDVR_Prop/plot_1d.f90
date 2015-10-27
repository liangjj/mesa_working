!deck plot_1d.f
!***begin prologue     plot_1d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_1d
  SUBROUTINE plot_1d(vector)
  USE arnoldi_global
  USE dvr_shared  
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: vector
  REAL*8                                 :: rtemp
  INTEGER                                :: i
  DO i = 1,nphy(1)
     rtemp = SQRT ( ( vector(i,1) * vector(i,1)        &
                                  +                      &
                      vector(i,2) * vector(i,2) ) )
     write(iplot(3),*) grid(1)%pt(i), vector(i,1)
     write(iplot(4),*) grid(1)%pt(i), vector(i,2)
     write(iplot(5),*) grid(1)%pt(i), rtemp
  END DO
END SUBROUTINE plot_1d
