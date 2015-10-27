!deck plot_psi_3d.f
!***begin prologue     plot_psi_3d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_psi_3d
  SUBROUTINE plot_psi_3d(vector,t)
  USE dvrprop_global_rt
  USE dvr_shared  
  USE dvr_global
  USE plot_wavefunction
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: vector
  REAL*8                                         :: rtemp
  INTEGER                                        :: i, j, k, t
  CHARACTER (LEN=16)                             :: fptoc
  IF (log_main(8)) THEN
      title='solution at t = '//fptoc(tim_pts(t+1))
      CALL print_psi(vector)
  END IF
  IF (plot.and.t == ntreg) THEN
      DO i = 1,nphy(3)
         DO j = 1,nphy(2)
            DO k = 1,nphy(1)
               rtemp = SQRT ( ( vector(i,j,k,1) * vector(i,j,k,1)      &
                                                +                      &
                                vector(i,j,k,2) * vector(i,j,k,2) ) )
               write(iplot(3),*) grid(3)%pt(i), grid(2)%pt(j),         &
                                 grid(1)%pt(k), vector(i,j,k,1)      
               write(iplot(4),*) grid(3)%pt(i), grid(2)%pt(j),         &
                                 grid(3)%pt(k), vector(i,j,k,2)
               write(iplot(5),*) grid(3)%pt(i), grid(2)%pt(j),         &
                                 grid(3)%pt(k), rtemp
            END DO
         END DO
      END DO
   END IF
END SUBROUTINE plot_psi_3d
