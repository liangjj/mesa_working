SUBROUTINE gauss_packet_1d(wave_function,norm)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)              :: wave_function
  REAL*8                                    :: norm
  INTEGER                                   :: i  
  ALLOCATE(zloc(1))
  ALLOCATE(zloc(1)%z(nphy(1),2))
  WRITE(iout,1)
  wave_function(:,:) = 0.d0
  IF(system == 'cartesian') THEN
     WRITE(iout,2)
     WRITE(iout,3)
     WRITE(iout,4) alpha(1)
     WRITE(iout,5) sigma(1)
     WRITE(iout,6) x_0(1)
     WRITE(iout,7) beta(1)
     CALL z_proj(grid(1)%pt,grid(1)%wt,grid(1)%f, &
                 alpha(1),sigma(1),x_0(1),beta(1),zloc(1)%z,nphy(1))
     wave_function(:,:) = zloc(1)%z(:,:)
  ELSE
     WRITE(iout,8)
     WRITE(iout,9) alpha(1), sigma(1), x_0(1), beta(1)
     CALL z_proj(grid(1)%pt,grid(1)%wt,grid(1)%f, &
                 alpha(1),sigma(1),x_0(1),beta(1),zloc(1)%z,nphy(1),  &
                 typke,system,log_main(4))
     wave_function(:,:) = zloc(1)%z(:,:)
  END IF
  wave_function(:,:) = norm * wave_function(:,:)
  DEALLOCATE(zloc(1)%z)
  DEALLOCATE(zloc)
1 FORMAT(/,5X,'initial wavepacket at t=0')
2 FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,  &
 'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
           ' - i*beta_x*(x-x_0) ]' )
3    FORMAT(/,1X,'gaussian wave packet parameters')
4    FORMAT(/,1X,'alpha    = ',e15.8)
5    FORMAT(/,1X,'sigma    = ',e15.8)
6    FORMAT(/,1X,'r_0      = ',e15.8)
7    FORMAT(/,1X,'beta     = ',e15.8)
8    FORMAT(/,5X,'the form of the radial packet is:',///,15X,  &
    'psi = r*exp(- alpha_r * ( r - r_0 )**2 / (2*(sigma_r)**2) - &
                               i*beta_r*(r-r_0))')
9    FORMAT(/,1X,'gaussian wave packet parameters'  &
           ,/,1X,'alpha = ',e15.8,/,1x,'sigma = ',e15.8, &
            /,1X,'r_0 = ',e15.8, /,1X,'beta  = ',e15.8)
END SUBROUTINE gauss_packet_1d
