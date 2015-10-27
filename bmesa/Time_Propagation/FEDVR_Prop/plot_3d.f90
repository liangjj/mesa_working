!deck plot_3d.f
!***begin prologue     plot_3d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           auto correlation function
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_3d
  SUBROUTINE plot_3d(vector)
  USE arnoldi_global
  USE dvr_shared  
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: vector
  REAL*8                                         :: rtemp
  INTEGER                                        :: i, j, k
  DO i = 1,nphy(3)
     DO j = 1,nphy(2)
        DO k = 1,nphy(1)
           rtemp = SQRT ( ( vector(i,j,k,1) * vector(i,j,k,1)      &
                                            +                      &
                            vector(i,j,k,2) * vector(i,j,k,2) ) )
          write(iplot(3),*) grid(3)%pt(i), grid(2)%pt(j),          &
                            grid(1)%pt(k), vector(i,j,k,1)      
          write(iplot(4),*) grid(3)%pt(i), grid(2)%pt(j),          &
                            grid(3)%pt(k), vector(i,j,k,2)
          write(iplot(5),*) grid(3)%pt(i), grid(2)%pt(j),          &
                            grid(3)%pt(k), rtemp
        END DO
     END DO
  END DO
END SUBROUTINE plot_3d
