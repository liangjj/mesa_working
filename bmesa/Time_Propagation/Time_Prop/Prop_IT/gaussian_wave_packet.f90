!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              MODULE gaussian_wave_packet
                              USE dvrprop_global
                              USE dvr_shared
                              USE dvr_global
                              USE z_project
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**begin prologue     gaussian_wave_packet
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
              MODULE PROCEDURE gauss_packet_d,                          &
                               gauss_packet_z
                          END INTERFACE gauss_packet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE gauss_packet_d(wave_function,norm)
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d)         :: wave_function
  REAL*8                         :: norm
  INTEGER                        :: i 
  ALLOCATE(zloc(spdim))
  DO i=1, spdim
     ALLOCATE(zloc(i)%z_d(nphy(i)))
  END DO
  WRITE(iout,1)
  wave_function(:) = 0.d0
  IF ( spdim == 1 ) THEN
       WRITE(iout,2)
       WRITE(iout,3)
       WRITE(iout,4) alpha(1)
       WRITE(iout,5) sigma(1)
       WRITE(iout,6) x_0(1)
       WRITE(iout,7) powr(1)
       WRITE(iout,8) beta(1)
  ELSE IF (spdim == 2) THEN
       WRITE(iout,9)
       WRITE(iout,3)
       WRITE(iout,10) (alpha(i),i=1,2)
       WRITE(iout,11) (sigma(i),i=1,2)
       WRITE(iout,12) (x_0(i),i=1,2)
       WRITE(iout,13) (powr(i),i=1,2)
       WRITE(iout,14) (beta(i),i=1,2)
  ELSE IF (spdim == 3) THEN
       WRITE(iout,15)
       WRITE(iout,3)
       WRITE(iout,16) (alpha(i),i=1,2)
       WRITE(iout,17) (sigma(i),i=1,2)
       WRITE(iout,18) (x_0(i),i=1,2)
       WRITE(iout,19) (powr(i),i=1,2)
       WRITE(iout,20) (beta(i),i=1,2)
  END IF
  DO i=1,spdim
     CALL z_proj(grid(i)%pt,grid(i)%wt,                                    &
                 alpha(i),sigma(i),x_0(i),powr(i),beta(i),                 &
                 zloc(i)%z_d,nphy(i))
  END DO
  IF (spdim == 1 ) THEN
      CALL wfn_fil_1d_dd(wave_function,zloc(1)%z_d)
  ELSE IF( spdim == 2) THEN
      CALL wfn_fil_2d_dd(wave_function,zloc(1)%z_d,zloc(2)%z_d)
  ELSE IF( spdim == 3) THEN
      CALL wfn_fil_3d_dd(wave_function,zloc(1)%z_d,zloc(2)%z_d,zloc(3)%z_d)
  END IF
  wave_function(:) = norm * wave_function(:)
  DO i=1,spdim
     DEALLOCATE(zloc(i)%z_d)
  END DO
  DEALLOCATE(zloc)
1  FORMAT(/,5X,'initial wavepacket at t=0')
2  FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,    &
    'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
              ' - i*beta_x*(x-x_0) ]' )
3 FORMAT(/,1X,'gaussian wave packet parameters')
4 FORMAT(/,1X,'alpha_x    = ',e15.8)
5 FORMAT(/,1X,'sigma_x    = ',e15.8)
6 FORMAT(/,1X,'x_0        = ',e15.8)
7 FORMAT(/,1X,'n_x        = ',i2)
8 FORMAT(/,1X,'beta_x     = ',e15.8)
9 FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,    &
    'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
              ' - i*beta_x*(x-x_0) ]',/,5X,                      &
    '      exp[ - alpha_y * ( y - y_0 )**2 / (2*(sigma_y)**2)'   &
              ' - i*beta_y*(y-y_0) ]' )
10 FORMAT(/,1X,'alpha_x    = ',e15.8,1X,'alpha_y   = ',e15.8)
11 FORMAT(/,1X,'sigma_x    = ',e15.8,1X,'sigma_y   = ',e15.8)
12 FORMAT(/,1X,'x_0        = ',e15.8,1X,'y_0       = ',e15.8)
13 FORMAT(/,1X,'n_x        = ',i2,   1X,'n_y       = ',i2 )
14 FORMAT(/,1X,'beta_x     = ',e15.8,1X,'beta_y    = ',e15.8 )
15 FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,    &
    'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
              ' - i*beta_x*(x-x_0) ]',/,5X,                      &
    '      exp[ - alpha_y * ( y - y_0 )**2 / (2*(sigma_y)**2)'   &
              ' - i*beta_y*(y-y_0) ]',/,5x,                      &
    '      exp[ - alpha_z * ( z - z_0 )**2 / (2*(sigma_z)**2)'   &
              ' - i*beta_z*(z-z_0) ]' )
16 FORMAT(/,1X,'alpha_x    = ',e15.8,1X,'alpha_y   = ',e15.8,1X,'alpha_z   = ',e15.8)
17 FORMAT(/,1X,'sigma_x    = ',e15.8,1X,'sigma_y   = ',e15.8,1X,'sigma_z   = ',e15.8)
18 FORMAT(/,1X,'x_0        = ',e15.8,1X,'y_0       = ',e15.8,1X,'z_0       = ',e15.8)
19 FORMAT(/,1X,'n_x        = ',i2,   1X,'n_y       = ',i2,   1X,'n_y       = ',i2)
20 FORMAT(/,1X,'beta_x     = ',e15.8,1X,'beta_y    = ',e15.8,'beta_z       = ',e15.8)
  END SUBROUTINE gauss_packet_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE gauss_packet_z(wave_function,norm)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d)         :: wave_function
  REAL*8                             :: norm
  INTEGER                            :: i  
  ALLOCATE(zloc(spdim))
  DO i=1, spdim
     ALLOCATE(zloc(i)%z_z(nphy(i)))
  END DO
  WRITE(iout,1)
  wave_function(:) = 0.d0
  IF ( spdim == 1 ) THEN
       WRITE(iout,2)
       WRITE(iout,3)
       WRITE(iout,4) alpha(1)
       WRITE(iout,5) sigma(1)
       WRITE(iout,6) x_0(1)
       WRITE(iout,7) powr(1)
       WRITE(iout,8) beta(1)
  ELSE IF (spdim == 2) THEN
       WRITE(iout,9)
       WRITE(iout,3)
       WRITE(iout,10) (alpha(i),i=1,2)
       WRITE(iout,11) (sigma(i),i=1,2)
       WRITE(iout,12) (x_0(i),i=1,2)
       WRITE(iout,13) (powr(i),i=1,2)
       WRITE(iout,14) (beta(i),i=1,2)
  ELSE IF (spdim == 3) THEN
       WRITE(iout,15)
       WRITE(iout,3)
       WRITE(iout,16) (alpha(i),i=1,2)
       WRITE(iout,17) (sigma(i),i=1,2)
       WRITE(iout,18) (x_0(i),i=1,2)
       WRITE(iout,19) (powr(i),i=1,2)
       WRITE(iout,20) (beta(i),i=1,2)
  END IF
  DO i=1,spdim
     CALL z_proj(grid(i)%pt,grid(i)%wt,                                    &
                 alpha(i),sigma(i),x_0(i),powr(i),beta(i),                 &
                 zloc(i)%z_z,nphy(i))
  END DO
  IF (spdim == 1 ) THEN
      CALL wfn_fil_1d_zz(wave_function,zloc(1)%z_z)
  ELSE IF( spdim == 2) THEN
      CALL wfn_fil_2d_zz(wave_function,zloc(1)%z_z,zloc(2)%z_z)
  ELSE IF( spdim == 3) THEN
      CALL wfn_fil_3d_zz(wave_function,zloc(1)%z_z,zloc(2)%z_z,zloc(3)%z_z)
  END IF
  wave_function(:) = norm * wave_function(:)
  DO i=1,spdim
     DEALLOCATE(zloc(i)%z_z)
  END DO
  DEALLOCATE(zloc)
1  FORMAT(/,5X,'initial wavepacket at t=0')
2  FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,    &
    'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
              ' - i*beta_x*(x-x_0) ]' )
3 FORMAT(/,1X,'gaussian wave packet parameters')
4 FORMAT(/,1X,'alpha_x    = ',e15.8)
5 FORMAT(/,1X,'sigma_x    = ',e15.8)
6 FORMAT(/,1X,'x_0        = ',e15.8)
7 FORMAT(/,1X,'n_x        = ',i2)
8 FORMAT(/,1X,'beta_x     = ',e15.8)
9 FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,    &
    'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
              ' - i*beta_x*(x-x_0) ]',/,5X,                      &
    '      exp[ - alpha_y * ( y - y_0 )**2 / (2*(sigma_y)**2)'   &
              ' - i*beta_y*(y-y_0) ]' )
10 FORMAT(/,1X,'alpha_x    = ',e15.8,1X,'alpha_y   = ',e15.8)
11 FORMAT(/,1X,'sigma_x    = ',e15.8,1X,'sigma_y   = ',e15.8)
12 FORMAT(/,1X,'x_0        = ',e15.8,1X,'y_0       = ',e15.8)
13 FORMAT(/,1X,'n_x        = ',i2,   1X,'n_y       = ',i2 )
14 FORMAT(/,1X,'beta_x     = ',e15.8,1X,'beta_y    = ',e15.8 )
15 FORMAT(/,5X,'the form of the cartesian packet is:',///,5X,    &
    'psi = exp[ - alpha_x * ( x - x_0 )**2 / (2*(sigma_x)**2)'   &
              ' - i*beta_x*(x-x_0) ]',/,5X,                      &
    '      exp[ - alpha_y * ( y - y_0 )**2 / (2*(sigma_y)**2)'   &
              ' - i*beta_y*(y-y_0) ]',/,5x,                      &
    '      exp[ - alpha_z * ( z - z_0 )**2 / (2*(sigma_z)**2)'   &
              ' - i*beta_z*(z-z_0) ]' )
16 FORMAT(/,1X,'alpha_x    = ',e15.8,1X,'alpha_y   = ',e15.8,1X,'alpha_z   = ',e15.8)
17 FORMAT(/,1X,'sigma_x    = ',e15.8,1X,'sigma_y   = ',e15.8,1X,'sigma_z   = ',e15.8)
18 FORMAT(/,1X,'x_0        = ',e15.8,1X,'y_0       = ',e15.8,1X,'z_0       = ',e15.8)
19 FORMAT(/,1X,'n_x        = ',i2,   1X,'n_y       = ',i2,   1X,'n_y       = ',i2)
20 FORMAT(/,1X,'beta_x     = ',e15.8,1X,'beta_y    = ',e15.8,'beta_z       = ',e15.8)
  END SUBROUTINE gauss_packet_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_1d_dd(wave_function,z_1)
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))            :: wave_function
  REAL*8, DIMENSION(nphy(1))            :: z_1
  wave_function(:) = z_1(:)
  END SUBROUTINE wfn_fil_1d_dd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_1d_zd(wave_function,z_1)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(1))        :: wave_function
  REAL*8, DIMENSION(nphy(1))            :: z_1
  wave_function(:) = z_1(:)
  END SUBROUTINE wfn_fil_1d_zd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_1d_zz(wave_function,z_1)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(1))            :: wave_function
  COMPLEX*16, DIMENSION(nphy(1))            :: z_1
  wave_function(:) = z_1(:)
  END SUBROUTINE wfn_fil_1d_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_2d_dd(wave_function,z_1,z_2)
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1))            :: wave_function
  REAL*8, DIMENSION(nphy(1))                    :: z_1
  REAL*8, DIMENSION(nphy(2))                    :: z_2
  INTEGER                                       :: i, j
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        wave_function(j,i) = z_1(i) * z_2(j)
     END DO
  END DO
  END SUBROUTINE wfn_fil_2d_dd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_2d_zd(wave_function,z_1,z_2)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(2),nphy(1))        :: wave_function
  REAL*8, DIMENSION(nphy(1))                    :: z_1
  REAL*8, DIMENSION(nphy(2))                    :: z_2
  INTEGER                                       :: i, j
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        wave_function(j,i) = z_1(i) * z_2(j)
     END DO
  END DO
  END SUBROUTINE wfn_fil_2d_zd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_2d_zz(wave_function,z_1,z_2)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(2),nphy(1))        :: wave_function
  COMPLEX*16, DIMENSION(nphy(1))                :: z_1
  COMPLEX*16, DIMENSION(nphy(2))                :: z_2
  INTEGER                                       :: i, j
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        wave_function(j,i) = z_1(i) * z_2(j)
     END DO
  END DO
  END SUBROUTINE wfn_fil_2d_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_3d_dd(wave_function,z_1,z_2,z_3)
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))    :: wave_function
  REAL*8, DIMENSION(nphy(1))                    :: z_1
  REAL*8, DIMENSION(nphy(2))                    :: z_2
  REAL*8, DIMENSION(nphy(3))                    :: z_3
  INTEGER                                       :: i, j, k
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        DO k=1,nphy(3)
           wave_function(k,j,i) = z_1(i) * z_2(j) * z_3(k)
        END DO
     END DO
  END DO
  END SUBROUTINE wfn_fil_3d_dd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_3d_zd(wave_function,z_1,z_2,z_3)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1))    :: wave_function
  REAL*8, DIMENSION(nphy(1))                        :: z_1
  REAL*8, DIMENSION(nphy(2))                        :: z_2
  REAL*8, DIMENSION(nphy(3))                        :: z_3
  INTEGER                                           :: i, j, k
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        DO k=1,nphy(3)
           wave_function(k,j,i) = z_1(i) * z_2(j) * z_3(k)
        END DO
     END DO
  END DO
  END SUBROUTINE wfn_fil_3d_zd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wfn_fil_3d_zz(wave_function,z_1,z_2,z_3)
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1))  :: wave_function
  COMPLEX*16, DIMENSION(nphy(1))                  :: z_1
  COMPLEX*16, DIMENSION(nphy(2))                  :: z_2
  COMPLEX*16, DIMENSION(nphy(3))                  :: z_3
  INTEGER                                       :: i, j, k
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        DO k=1, nphy(3)
           wave_function(k,j,i) = z_1(i) * z_2(j) * z_3(k)
        END DO
     END DO
  END DO
  END SUBROUTINE wfn_fil_3d_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END MODULE gaussian_wave_packet
