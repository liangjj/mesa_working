! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Subroutine for Gaussian Wavepacket}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck gauss_paket.f
!**begin prologue     gauss_paket
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate zero time gaussian wavepacket
!**references
!**routines called
!**end prologue       gauss_paket
  SUBROUTINE gauss_paket(norm)
  USE arnoldi_global
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  REAL*8                                 :: norm
  COMPLEX*16, DIMENSION(:), ALLOCATABLE  :: z
  INTEGER                                :: i  
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
! Go get some needed scratch memory.
  ALLOCATE(z(n3d))
  WRITE(iout,1)
  psi0 = (0.d0,0.d0)
  IF(system == 'cartesian') THEN
     WRITE(iout,2)
     WRITE(iout,3)
     WRITE(iout,4) (alpha(i),i=1,spdim)
     WRITE(iout,5) (sigma(i),i=1,spdim)
     WRITE(iout,6) (x_0(i),i=1,spdim)
     WRITE(iout,7) (beta(i),i=1,spdim)
     IF(spdim == 1) THEN
        CALL z_proj(grid(1)%pt,grid(1)%wt,grid(1)%f, &
                    alpha(1),sigma(1),x_0(1),beta(1),z,nphy(1))
        psi0 = z(1:n3d)
     ELSE IF(spdim == 2) THEN
        CALL z_proj(grid(1)%pt,grid(1)%wt,grid(1)%f, &
                    alpha(1),sigma(1),x_0(1),beta(1),z,nphy(1))
        CALL z_proj(grid(2)%pt,grid(2)%wt,grid(2)%f, &
                    alpha(2),sigma(2),x_0(2),beta(2),z(1+nphy(1)), &
                    nphy(2))
        CALL zfil_2(z,z(1+nphy(1)))
     ELSE IF(spdim == 3) THEN
        CALL z_proj(grid(1)%pt,grid(1)%wt,grid(1)%f, &
                    alpha(1),sigma(1),x_0(1),beta(1),z,nphy(1))
        CALL z_proj(grid(2)%pt,grid(2)%wt,grid(2)%f, &
                    alpha(2),sigma(2),x_0(2),beta(2),z(1+nphy(1)), &
                    nphy(2))    
        CALL z_proj(grid(3)%pt,grid(3)%wt,grid(3)%f, &
                    alpha(3),sigma(3),x_0(3),beta(3),z(1+nphy(1)+nphy(2)), &
                    nphy(3))
        CALL zfil_3(z,z(1+nphy(1)),z(1+nphy(1)+nphy(2)))
    END IF
  ELSE
    WRITE(iout,8)
    WRITE(iout,9) alpha(1), sigma(1), x_0(1), beta(1)
    CALL z_proj(grid(1)%pt,grid(1)%wt,grid(1)%f, &
                alpha(1),sigma(1),x_0(1),beta(1),z,nphy(1),  &
                typke,system,log_main(4))
    psi0 = z(1:n3d)
  END IF
  psi0 = norm*psi0
  DEALLOCATE(z)
1    FORMAT(/,5X,'initial wavepacket at t=0')
2    FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,  &
    'psi = exp( - alpha_x * ( x - x_0 )**2 - i*beta_x*(x-x_0) )  &
                                      / (2*(sigma_x)**2)',/,5X,  &
    '      exp( - alpha_y * ( y - y_0 )**2 - i*beta_y*(y-y_0) )  &
                                      / (2*(sigma_y)**2)',/,5X,  &
    '      exp( - alpha_z * ( z - z_0 )**2 - i*beta_z*(z-z_0) )  &
                                      / (2*(sigma_z)**2)')
3    FORMAT(/,1X,'gaussian wave packet parameters')
4    FORMAT(/,1X,'alpha    = ',3(e15.8,1X))
5    FORMAT(/,1X,'sigma    = ',3(e15.8,1X))
6    FORMAT(/,1X,'r_0      = ',3(e15.8,1X))
7    FORMAT(/,1X,'beta     = ',3(e15.8,1X))
8    FORMAT(/,5X,'the form of the radial packet is:',///,15X,  &
    'psi = r*exp(- alpha_r * ( r - r_0 )**2 / (2*(sigma_r)**2) - &
                               i*beta_r*(r-r_0))')
9    FORMAT(/,1X,'gaussian wave packet parameters'  &
           ,/,1X,'alpha = ',e15.8,/,1x,'sigma = ',e15.8, &
            /,1X,'r_0 = ',e15.8, /,1X,'beta  = ',e15.8)
END SUBROUTINE gauss_paket







