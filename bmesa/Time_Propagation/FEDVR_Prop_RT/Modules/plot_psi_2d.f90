!deck plot_psi_2d.f
!***begin prologue     plot_psi_2d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_psi_2d
  SUBROUTINE plot_psi_2d(vector,t)
  USE dvrprop_global_rt
  USE dvr_shared  
  USE dvr_global
  USE plot_wavefunction
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: vector
  REAL*8                                 :: rtemp
  INTEGER                                :: t, i, j
  CHARACTER (LEN=16)                     :: fptoc
  IF (log_main(8)) THEN
      title='solution at t = '//fptoc(tim_pts(t+1))
      CALL print_psi(vector)
  END IF
  IF (plot.and.t == ntreg) THEN
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
  END IF
END SUBROUTINE plot_psi_2d
