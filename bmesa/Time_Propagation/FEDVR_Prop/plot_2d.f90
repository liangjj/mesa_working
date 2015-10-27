!deck plot_2d.f
!***begin prologue     plot_2d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_2d
  SUBROUTINE plot_2d(vector)
  USE arnoldi_global
  USE dvr_shared  
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: vector
  REAL*8                                 :: rtemp
  INTEGER                                :: i, j
  DO i = 1,nphy(2)
     DO j = 1,nphy(1)
        rtemp = SQRT ( ( vector(i,j,1) * vector(i,j,1)        &
                                       +                      &
                         vector(i,j,2) * vector(i,j,2) ) )
        write(iplot(3),*) grid(2)%pt(i), grid(1)%pt(j), vector(i,j,1)
        write(iplot(4),*) grid(2)%pt(i), grid(1)%pt(j), vector(i,j,2)
        write(iplot(5),*) grid(2)%pt(i), grid(1)%pt(j), rtemp
     END DO
  END DO
END SUBROUTINE plot_2d
