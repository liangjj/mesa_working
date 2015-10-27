!deck plot_psi_1d.f
!***begin prologue     plot_psi_1d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose             
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       plot_psi_1d
  SUBROUTINE plot_psi_1d(vector,t)
  USE dvrprop_global_rt
  USE dvr_shared  
  USE dvr_global
  USE plot_wavefunction
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: vector
  REAL*8                                 :: rtemp
  INTEGER                                :: t, i
  CHARACTER (LEN=16)                     :: fptoc
  IF (log_main(8)) THEN
      title='solution at t = '//fptoc(tim_pts(t+1))
      CALL print_psi(vector)
  END IF
  IF (plot.and.t == ntreg) THEN
      DO i = 1,nphy(1)
         rtemp = SQRT ( ( vector(i,1) * vector(i,1)        &
                                      +                    &
                          vector(i,2) * vector(i,2) ) )
         write(iplot(3),*) grid(1)%pt(i), vector(i,1)
         write(iplot(4),*) grid(1)%pt(i), vector(i,2)
         write(iplot(5),*) grid(1)%pt(i), rtemp
      END DO
  END IF
END SUBROUTINE plot_psi_1d
