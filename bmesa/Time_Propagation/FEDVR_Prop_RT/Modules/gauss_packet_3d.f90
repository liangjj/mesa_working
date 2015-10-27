SUBROUTINE gauss_packet_3d(wave_function,norm)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)       :: wave_function
  REAL*8                                             :: norm
  INTEGER                                            :: i, j, k  
  ALLOCATE(zloc(3))
  ALLOCATE(zloc(1)%z(nphy(1),2),zloc(2)%z(nphy(2),2),zloc(3)%z(nphy(3),2))
  WRITE(iout,1)
  wave_function(:,:,:,:) = 0.d0
  WRITE(iout,2)
  WRITE(iout,3)
  WRITE(iout,4) (alpha(i),i=1,3)
  WRITE(iout,5) (sigma(i),i=1,3)
  WRITE(iout,6) (x_0(i),i=1,3)
  WRITE(iout,7) (beta(i),i=1,3)
  DO i=1,3
     CALL z_proj(grid(i)%pt,grid(i)%wt,grid(i)%f,                                    &
                 alpha(i),sigma(i),x_0(i),beta(i),zloc(i)%z,nphy(i))
  END DO
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        DO k=1,nphy(3)
           wave_function(k,j,i,1) = zloc(1)%z(i,1) * zloc(2)%z(j,1) * zloc(3)%z(k,1) &
                                                   -                                 &
                                    zloc(1)%z(i,2) * zloc(2)%z(j,2) * zloc(3)%z(k,1) &
                                                   -                                 &
                                    zloc(1)%z(i,2) * zloc(2)%z(j,1) * zloc(3)%z(k,2) &
                                                   -                                 &
                                    zloc(1)%z(i,1) * zloc(2)%z(j,2) * zloc(3)%z(k,2)
           wave_function(k,j,i,2) = zloc(1)%z(i,2) * zloc(2)%z(j,1) * zloc(3)%z(k,1) &
                                                   +                                 &
                                    zloc(1)%z(i,1) * zloc(2)%z(j,2) * zloc(3)%z(k,1) &
                                                   +                                 &
                                    zloc(1)%z(i,1) * zloc(2)%z(j,1) * zloc(3)%z(k,2) &
                                                   -                                 &
                                    zloc(1)%z(i,2) * zloc(2)%z(j,2) * zloc(3)%z(k,2)
        END DO
     END DO
  END DO
  wave_function(:,:,:,:) = norm * wave_function(:,:,:,:)
  DEALLOCATE(zloc(1)%z,zloc(2)%z,zloc(3)%z)
  DEALLOCATE(zloc)
1    FORMAT(/,5X,'initial wavepacket at t=0')
2    FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,  &
    'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
              ' - i*beta_x*(x-x_0) ]',/,5X,                      &
    '      exp[ - alpha_y * ( y - y_0 )**2 / (2*(sigma_y)**2)'   &
              ' - i*beta_y*(y-y_0) ]',/,5x,                      &
    '      exp[ - alpha_z * ( z - z_0 )**2 / (2*(sigma_z)**2)'   &
              ' - i*beta_z*(z-z_0) ]' )
3    FORMAT(/,1X,'gaussian wave packet parameters')
4    FORMAT(/,1X,'alpha    = ',3(e15.8,1X))
5    FORMAT(/,1X,'sigma    = ',3(e15.8,1X))
6    FORMAT(/,1X,'r_0      = ',3(e15.8,1X))
7    FORMAT(/,1X,'beta     = ',3(e15.8,1X))
END SUBROUTINE gauss_packet_3d
