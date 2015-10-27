!deck rad_paket_08
!**begin prologue     rad_paket_08
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate zero time rad wavepacket
!**references
!**routines called
!**end prologue       rad_paket_08
  SUBROUTINE rad_paket_08(psi)
  USE dvrprop_global
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d,2)               :: psi
  REAL*8                                 :: norm, sdot
  INTEGER                                :: i  
  WRITE(iout,1)
  IF(typ_pak == 'exponential') then
     WRITE(iout,2)
     psi(:,1) = ( grid(1)%pt - x_0(1) )**powr                    &
                                    *                            &
                exp( - ( sigma(1) ) * ( grid(1)%pt - x_0(1) ) )
     psi(:,2) = - psi(:,1) * sin( beta(1)  * ( grid(1)%pt - x_0(1) ) )
     psi(:,1) =   psi(:,1) * cos( beta(1)  * ( grid(1)%pt - x_0(1) ) )
  ELSE IF(typ_pak == 'gaussian') then
     WRITE(iout,3)
     psi(:,1) = ( grid(1)%pt - x_0(1) )**powr                    &
                                 *                               &
                exp( - sigma(1)  * ( grid(1)%pt - x_0(1) )       &
                                 * ( grid(1)%pt - x_0(1) ) )
     psi(:,2) = - sin ( grid(1)%pt - x_0(1) ) * psi(:,1)
     psi(:,1) =   cos ( grid(1)%pt - x_0(1) ) * psi(:,1)
  END IF
  WRITE(iout,4)
  WRITE(iout,5) sigma(1), x_0(1), beta(1), powr
  psi(:,1)  = psi(:,1) * grid(1)%wt
  psi(:,2)  = psi(:,2) * grid(1)%wt
  do i=1,nphy(1)
     psi(i,1) = psi(i,1) * grid(1)%f(i,i)
  END DO
  norm = 1.d0/sqrt( sdot(nphy(1),psi(1,1),1,psi(1,1),1)          &
                                  +                              &
                    sdot(nphy(1),psi(1,2),1,psi(1,2),1) )
  psi = psi * norm  
  norm = sdot(nphy(1),psi(1,1),1,psi(1,1),1)                     &
                                  +                              &
         sdot(nphy(1),psi(1,2),1,psi(1,2),1)
  title='initial wavepacket'
  CALL plot_wavefunction(psi)
  write(iout,7) norm
1    FORMAT(/,5X,'initial wavepacket at t=0')
2    FORMAT(/,5X,'the form of the radial packet is:',///,5X,     &
                 'psi = (r - x_0 )**n * '                        &
                 'exp( - ( alpha + i*beta ) * ( r - x_0 ))')
3    FORMAT(/,5X,'the form of the radial packet is:',///,5X,     &
                 'psi = ( r - x_0 )**n * '                       &
                 'exp( - alpha * ( r -x_0 ) * ( r - x_0 ) - '    &
                 'i * beta * (r- x_0 ) )')
4    FORMAT(/,1X,'radial wave packet parameters')
5    FORMAT(/,1X,'exponent    = ',e15.8,                         &
            /,1X,'shift       = ',e15.8,                         &
            /,1X,'momentum    = ',e15.8,                         &
            /,1X,'power       = ',i3)
7    FORMAT(/,2x,'Normalization of Initial Wavepacket = ',e15.8,1x,e15.8)
END SUBROUTINE rad_paket_08
