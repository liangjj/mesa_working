!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              MODULE gaussian_wave_packet
!**begin prologue     gaussian_wave_packett
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate zero time gaussian wavepacket
!**references
!**routines called
!**end prologue       gaussian_wave_packet
! \center Gaussian Wave Packet
! \par
!\begin{eqnarray}
! \Psi({\bf r}) & = N *
!                  &exp( - {\alpha}_x \frac{ ( x -x_0 )^2 }{2 { {\sigma}_x }^2 }
!                       - i {\beta}_x ( x - x_0 ) ) & \nonumber \\
!                  &exp( - {\alpha}_y \frac{ ( y -y_0 )^2 }{2 { {\sigma}_y }^2 }
!                        - i {\beta}_y ( y - y_0 ) ) & \nonumber \\
!                  &exp( - {\alpha}_z \frac{ ( z -z_0 )^2 }{2 { {\sigma}_z }^2 }
!                        - i {\beta}_z ( z - z_0 ) ) &
! \end{eqnarray}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              INTERFACE gauss_packet
              MODULE PROCEDURE gauss_packet_1d_d, gauss_packet_2d_d,   &
                               gauss_packet_3d_d
                          END INTERFACE gauss_packet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gauss_packet_1d_d(wave_function,norm)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))                :: wave_function
  REAL*8                                    :: norm
  INTEGER                                   :: i  
  ALLOCATE(zloc(1))
  ALLOCATE(zloc(1)%z(nphy(1)))
  WRITE(iout,1)
  wave_function = 0.d0
  IF(system == 'cartesian') THEN
     WRITE(iout,2)
     WRITE(iout,3)
     WRITE(iout,4) alpha(1)
     WRITE(iout,5) sigma(1)
     WRITE(iout,6) x_0(1)
     CALL z_proj(grid(1)%pt,grid(1)%wt,                       &
                 alpha(1),sigma(1),x_0(1),zloc(1)%z,nphy(1))
     wave_function = zloc(1)%z
  ELSE
     WRITE(iout,8)
     WRITE(iout,9) alpha(1), sigma(1), x_0(1)
     CALL z_proj(grid(1)%pt,grid(1)%wt,                       &
                 alpha(1),sigma(1),x_0(1),zloc(1)%z,nphy(1),  &
                 typke,system,log_main(4))
     wave_function = zloc(1)%z
  END IF
  wave_function = norm * wave_function
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
8    FORMAT(/,5X,'the form of the radial packet is:',///,15X,  &
    'psi = r*exp(- alpha_r * ( r - r_0 )**2 / (2*(sigma_r)**2) - &
                               i*beta_r*(r-r_0))')
9    FORMAT(/,1X,'gaussian wave packet parameters'  &
           ,/,1X,'alpha = ',e15.8,/,1x,'sigma = ',e15.8, &
            /,1X,'r_0 = ',e15.8)
END SUBROUTINE gauss_packet_1d_d
SUBROUTINE gauss_packet_2d_d(wave_function,norm)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1))        :: wave_function
  REAL*8                                    :: norm
  INTEGER                                   :: i, j  
  ALLOCATE(zloc(2))
  ALLOCATE(zloc(1)%z(nphy(1)),zloc(2)%z(nphy(2)))
  WRITE(iout,1)
  wave_function = 0.d0
  WRITE(iout,2)
  WRITE(iout,3)
  WRITE(iout,4) (alpha(i),i=1,2)
  WRITE(iout,5) (sigma(i),i=1,2)
  WRITE(iout,6) (x_0(i),i=1,2)
  DO i=1,2
     CALL z_proj(grid(i)%pt,grid(i)%wt,                      &
                 alpha(i),sigma(i),x_0(i),zloc(i)%z,         &
                 nphy(i))
  END DO
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        wave_function(j,i) = zloc(1)%z(i) * zloc(2)%z(j)   
     END DO
  END DO
  wave_function = norm * wave_function
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
END SUBROUTINE gauss_packet_2d_d
SUBROUTINE gauss_packet_3d_d(wave_function,norm)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))         :: wave_function
  REAL*8                                             :: norm
  INTEGER                                            :: i, j, k  
  ALLOCATE(zloc(3))
  ALLOCATE(zloc(1)%z(nphy(1)),zloc(2)%z(nphy(2)),zloc(3)%z(nphy(3)))
  WRITE(iout,1)
  wave_function = 0.d0
  WRITE(iout,2)
  WRITE(iout,3)
  WRITE(iout,4) (alpha(i),i=1,3)
  WRITE(iout,5) (sigma(i),i=1,3)
  WRITE(iout,6) (x_0(i),i=1,3)
  DO i=1,3
     CALL z_proj(grid(i)%pt,grid(i)%wt,                           &
                 alpha(i),sigma(i),x_0(i),zloc(i)%z,nphy(i))
  END DO
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        DO k=1,nphy(3)
           wave_function(k,j,i) = zloc(1)%z(i)                    &
                                     *                            &
                                  zloc(2)%z(j)                    &
                                     *                            &
                                  zloc(3)%z(k)
        END DO
     END DO
  END DO
  wave_function = norm * wave_function
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
END SUBROUTINE gauss_packet_3d_d
END MODULE gaussian_wave_packet
