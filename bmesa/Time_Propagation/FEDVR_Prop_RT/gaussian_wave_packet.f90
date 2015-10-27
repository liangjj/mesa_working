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
             MODULE PROCEDURE gauss_packet_1d_z, gauss_packet_2d_z,    &
                              gauss_packet_3d_z
                     END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gauss_packet_1d_z(wave_function,norm)
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
     CALL z_proj(grid(1)%pt,grid(1)%wt,                              &
                 alpha(1),sigma(1),x_0(1),beta(1),zloc(1)%z,nphy(1))
     wave_function(:,:) = zloc(1)%z(:,:)
  ELSE
     WRITE(iout,8)
     WRITE(iout,9) alpha(1), sigma(1), x_0(1), beta(1)
     CALL z_proj(grid(1)%pt,grid(1)%wt,                               &
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
6    FORMAT(/,1X,'x_0      = ',e15.8)
7    FORMAT(/,1X,'beta     = ',e15.8)
8    FORMAT(/,5X,'the form of the radial packet is:',///,15X,  &
    'psi = r*exp(- alpha_r * ( r - r_0 )**2 / (2*(sigma_r)**2) - &
                               i*beta_r*(r-r_0))')
9    FORMAT(/,1X,'gaussian wave packet parameters'  &
           ,/,1X,'alpha = ',e15.8,/,1x,'sigma = ',e15.8, &
            /,1X,'r_0 = ',e15.8, /,1X,'beta  = ',e15.8)
END SUBROUTINE gauss_packet_1d_z
SUBROUTINE gauss_packet_2d_z(wave_function,norm)
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
     CALL z_proj(grid(i)%pt,grid(i)%wt,                           &
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
6 FORMAT(/,1X,'x_0      = ',2(e15.8,1X))
7 FORMAT(/,1X,'beta     = ',2(e15.8,1X))
END SUBROUTINE gauss_packet_2d_z
SUBROUTINE gauss_packet_3d_z(wave_function,norm)
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
     CALL z_proj(grid(i)%pt,grid(i)%wt,                                &
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
6    FORMAT(/,1X,'x_0      = ',3(e15.8,1X))
7    FORMAT(/,1X,'beta     = ',3(e15.8,1X))
END SUBROUTINE gauss_packet_3d_z
! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Subroutine for Gaussian Wavepacket Normalization}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck nr_packet.f
!**begin prologue     nr_packet
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate normailzation integrals for gaussian wavepacket
!**references
!**routines called
!**end prologue       nr_packet
  SUBROUTINE nr_packet(norm)
  USE arnoldi_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8                                 :: norm
  REAL*8, DIMENSION(3)                   :: ov 
  INTEGER                                :: i
  IF(system == 'cartesian') THEN
     norm=1.d0
     DO i=1,spdim
        call ov1_quad(ov(i),grid(i)%pt,grid(i)%wt, &
                      x_0(i),alpha(i),sigma(i),nphy(i))
        norm=norm*ov(i)
     END DO
  ELSE
     call ov2_quad(norm,grid(1)%pt,grid(1)%wt, &
                   x_0(1),alpha(1),sigma(1),nphy(1))    
  END IF
  norm=1.d0/SQRT(norm)
  WRITE(iout,1) norm
1    FORMAT(/,1X,'overlap integral for gaussian packet = ',e15.8)
END SUBROUTINE nr_packet
!deck rad_packet
!**begin prologue     rad_packet
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate zero time rad wavepacket
!**references
!**routines called
!**end prologue       rad_packet
  SUBROUTINE rad_packet(wave_function)
  USE dvrprop_global_rt
  USE dvr_global
  USE dvr_shared
  USE plot_wavefunction
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: wave_function
  REAL*8                                 :: norm, sdot
  INTEGER                                :: i  
  WRITE(iout,1)
  IF(typ_pak == 'exponential') then
     WRITE(iout,2)
     wave_function(:,1) = ( grid(1)%pt - x_0(1) )**powr                        &
                                       *                                       &
                       exp( - ( sigma(1) ) * ( grid(1)%pt - x_0(1) ) )
     wave_function(:,2) = - wave_function(:,1) * sin( beta(1)                  &
                                               *                               &
                                  ( grid(1)%pt - x_0(1) ) )
     wave_function(:,1) =   wave_function(:,1) * cos( beta(1)                  &
                                               * ( grid(1)%pt - x_0(1) ) )
  ELSE IF(typ_pak == 'gaussian') then
     WRITE(iout,3)
     wave_function(:,1) = ( grid(1)%pt - x_0(1) )**powr                        &
                                       *                                       &
                       exp( - sigma(1) * ( grid(1)%pt - x_0(1) )               &
                                       *                                       &
                             ( grid(1)%pt - x_0(1) ) )
                   
     wave_function(:,2) = - sin ( grid(1)%pt - x_0(1) ) * wave_function(:,1)
     wave_function(:,1) =   cos ( grid(1)%pt - x_0(1) ) * wave_function(:,1)
  END IF
  WRITE(iout,4)
  WRITE(iout,5) sigma(1), x_0(1), beta(1), powr
  wave_function(:,1)  = wave_function(:,1) * sqrt( grid(1)%wt )
  wave_function(:,2)  = wave_function(:,2) * sqrt (grid(1)%wt )
  norm = 1.d0/sqrt( sdot(nphy(1),wave_function(1,1),1,wave_function(1,1),1)    &
                                              +                                &
                    sdot(nphy(1),wave_function(1,2),1,wave_function(1,2),1) )
  wave_function = wave_function * norm  
  norm = sdot(nphy(1),wave_function(1,1),1,wave_function(1,1),1)               &
                                              +                                &
         sdot(nphy(1),wave_function(1,2),1,wave_function(1,2),1)
  title='initial wavepacket'
  CALL print_psi(wave_function)
  write(iout,7) norm
1    FORMAT(/,5X,'initial wavepacket at t=0')
2    FORMAT(/,5X,'the form of the radial packet is:',///,5X,  &
                 'psi = (r - x_0 )**n * exp( - ( alpha + i*beta ) * ' &
                                           '( r - x_0 ))')
3    FORMAT(/,5X,'the form of the radial packet is:',///,5X,  &
                 'psi = ( r - x_0 )**n * &
                  exp( - alpha * ( r -x_0 ) * ( r - x_0 ) - i * beta * ' &
                                                            '(r- x_0 ) )')
4    FORMAT(/,1X,'radial wave packet parameters')
5    FORMAT(/,1X,'exponent    = ',e15.8, &
            /,1X,'shift       = ',e15.8, &
            /,1X,'momentum    = ',e15.8, &
            /,1X,'power       = ',i3)
7    FORMAT(/,2x,'Normalization of Initial Wavepacket = ',e15.8,1x,e15.8)
END SUBROUTINE rad_packet
END MODULE gaussian_wave_packet
