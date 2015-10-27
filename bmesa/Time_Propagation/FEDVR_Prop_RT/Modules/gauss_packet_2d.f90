SUBROUTINE gauss_packet_2d(wave_function,norm)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),2)      :: wave_function
  REAL*8                                    :: norm
  INTEGER                                   :: i, j  
  ALLOCATE(zloc(2))
  ALLOCATE(zloc(1)%z(nphy(1),2),zloc(2)%z(nphy(2),2))
  WRITE(iout,1)
  wave_function = 0.d0
  WRITE(iout,2)
  WRITE(iout,3)
  WRITE(iout,4) (alpha(i),i=1,2)
  WRITE(iout,5) (sigma(i),i=1,2)
  WRITE(iout,6) (x_0(i),i=1,2)
  WRITE(iout,7) (beta(i),i=1,2)
  DO i=1,2
     CALL z_proj(grid(i)%pt,grid(i)%wt,grid(i)%f,                 &
                 alpha(i),sigma(i),x_0(i),beta(i),zloc(i)%z,      &
                 nphy(i))
  END DO
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        wave_function(j,i,1) = zloc(1)%z(i,1) * zloc(2)%z(j,1)    &
                                              -                   &
                               zloc(1)%z(i,2) * zloc(2)%z(j,2) 
        wave_function(j,i,2) = zloc(1)%z(i,1) * zloc(2)%z(j,2)    &
                                              +                   &
                               zloc(1)%z(i,2) * zloc(2)%z(j,1) 
     END DO
  END DO
  wave_function(:,:,:) = norm * wave_function(:,:,:)
  DEALLOCATE(zloc(1)%z,zloc(2)%z)
  DEALLOCATE(zloc)
1 FORMAT(/,5X,'initial wavepacket at t=0')
2 FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,  &
  'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'  &
            ' - i*beta_x*(x-x_0) ]',/,5X,                     &
 '       exp[ - alpha_y * ( y - y_0 )**2 / (2*(sigma_y)**2)'  &
            ' - i*beta_y*(y-y_0) ]' )
3 FORMAT(/,1X,'gaussian wave packet parameters')
4 FORMAT(/,1X,'alpha    = ',2(e15.8,1X))
5 FORMAT(/,1X,'sigma    = ',2(e15.8,1X))
6 FORMAT(/,1X,'r_0      = ',2(e15.8,1X))
7 FORMAT(/,1X,'beta     = ',2(e15.8,1X))
END SUBROUTINE gauss_packet_2d
