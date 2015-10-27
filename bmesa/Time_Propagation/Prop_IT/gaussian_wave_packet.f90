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
!               & { ( x -x _0 ) }^{n}_x exp( - {\alpha}_x \frac{ ( x -x_0 )^2 }
!                                                           {2 { {\sigma}_x }^2 } )
!                                                            & \nonumber \\
!               & { ( y - y_0 ) }^{n}_y exp( - {\alpha}_y \frac{ ( y -y_0 )^2 } 
!                                                           {2 { {\sigma}_y }^2 } )
!                                                            & \nonumber \\
!               & { (z -z _0 ) }^{n}_z exp( - {\alpha}_z \frac{ ( z -z_0 )^2 }
!                                                           {2 { {\sigma}_z }^2 } )
!                                                            &
! \end{eqnarray}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              INTERFACE gauss_packet
              MODULE PROCEDURE gauss_packet_1d, gauss_packet_2d,   &
                               gauss_packet_3d
                          END INTERFACE gauss_packet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gauss_packet_1d(wave_function,norm)
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
  WRITE(iout,2)
  WRITE(iout,3)
  WRITE(iout,4) alpha(1)
  WRITE(iout,5) sigma(1)
  WRITE(iout,6) x_0(1)
  WRITE(iout,7) powr(1)
  CALL z_proj(grid(1)%pt,grid(1)%wt,                                   &
              alpha(1),sigma(1),x_0(1),powr(1),zloc(1)%z,nphy(1))
  wave_function = norm * zloc(1)%z
  DEALLOCATE(zloc(1)%z)
  DEALLOCATE(zloc)
1 FORMAT(/,5X,'initial wavepacket at t=0')
2 FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,                &
              'psi =  ( x - x_0) ** n_x * exp[ - alpha_x * ( x - x_0 )**2 ' &
              '/ (2*(sigma_x)**2) ]' )
3    FORMAT(/,1X,'gaussian wave packet parameters')
4    FORMAT(/,1X,'alpha_x    = ',e15.8)
5    FORMAT(/,1X,'sigma_x    = ',e15.8)
6    FORMAT(/,1X,'x_0        = ',e15.8)
7    FORMAT(/,1X,'n_x        = ',i2)
END SUBROUTINE gauss_packet_1d
SUBROUTINE gauss_packet_2d(wave_function,norm)
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
  WRITE(iout,7) (powr(i),i=1,2)
  DO i=1,2
     CALL z_proj(grid(i)%pt,grid(i)%wt,                                &
                 alpha(i),sigma(i),x_0(i),powr(i),zloc(i)%z,nphy(i))
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
2 FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,                &
              'psi = ( x - x_0 ) ** n_x * exp[ - alpha_x * ( x - x_0 )**2 ' &
              '/ (2*(sigma_x)**2) ]'                    ,/,5x               &
              '      ( y - y_0 ) ** n_y * exp[ - alpha_y * ( y - y_0 )**2 ' &
              '/ (2*(sigma_y)**2) ]' )
3 FORMAT(/,1X,'gaussian wave packet parameters')
4 FORMAT(/,1X,'alpha_x    = ',e15.8,1X,'alpha_y   = ',e15.8)
5 FORMAT(/,1X,'sigma_x    = ',e15.8,1X,'sigma_y   = ',e15.8)
6 FORMAT(/,1X,'x_0        = ',e15.8,1X,'y_0       = ',e15.8)
7 FORMAT(/,1X,'n_x        = ',i2,   1X,'n_y       = ',i2)
END SUBROUTINE gauss_packet_2d
SUBROUTINE gauss_packet_3d(wave_function,norm)
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
  WRITE(iout,7) (powr(i),i=1,3)
  DO i=1,3
     CALL z_proj(grid(i)%pt,grid(i)%wt,                                &
                 alpha(i),sigma(i),x_0(i),powr(i),zloc(i)%z,nphy(i))
  END DO
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        DO k=1,nphy(3)
           wave_function(k,j,i) = zloc(1)%z(i)                         &
                                     *                                 &
                                  zloc(2)%z(j)                         &
                                     *                                 &
                                  zloc(3)%z(k)
        END DO
     END DO
  END DO
  wave_function = norm * wave_function
  DEALLOCATE(zloc(1)%z,zloc(2)%z,zloc(3)%z)
  DEALLOCATE(zloc)
1    FORMAT(/,5X,'initial wavepacket at t=0')
2    FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,        &
                 'psi = ( x - x_0 ) ** n_x * exp[ - alpha_x * ( x - x_0 )**2 ' &
                 '/ (2*(sigma_x)**2) ]'                    ,/,5x       &
                 '      ( y - y_0 ) ** n_y * exp[ - alpha_y * ( y - y_0 )**2 ' &
                 '/ (2*(sigma_y)**2) ]'                    ,/,5x               &
                 '      ( z - z_0 ) ** n_z * exp[ - alpha_z * ( z - z_0 )**2 ' &
                 '/ (2*(sigma_z)**2) ]' )
3 FORMAT(/,1X,'gaussian wave packet parameters')
4 FORMAT(/,1X,'alpha_x    = ',e15.8,1X,'alpha_y   = ',e15.8,1X,'alpha_z   = ',e15.8)
5 FORMAT(/,1X,'sigma_x    = ',e15.8,1X,'sigma_y   = ',e15.8,1X,'sigma_z   = ',e15.8)
6 FORMAT(/,1X,'x_0        = ',e15.8,1X,'y_0       = ',e15.8,1X,'z_0       = ',e15.8)
7 FORMAT(/,1X,'n_x        = ',i2,   1X,'n_y       = ',i2,   1X,'n_y       = ',i2)
END SUBROUTINE gauss_packet_3d
END MODULE gaussian_wave_packet
