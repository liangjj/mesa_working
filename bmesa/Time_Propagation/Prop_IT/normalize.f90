!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             MODULE normalize
!**begin prologue     normalize
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       normalize
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             INTERFACE check_norm
             MODULE PROCEDURE check_norm_1d_d, check_norm_2d_d,      &
                              check_norm_3d_d 
                             END INTERFACE check_norm
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck check_norm_1d_d.f
!***begin prologue     check_norm_1d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for one dimension wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_norm_1d_d
  SUBROUTINE check_norm_1d_d(v,norm)
  USE dvrprop_global_it
  USE dvr_shared,           ONLY  : nphy
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))             :: v
  REAL*8                                 :: ddot
  REAL*8                                 :: norm
!
  norm = ddot(nphy(1),v,1,v,1)
  write(iout,1) norm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE check_norm_1d_d
!deck check_norm_2d_d.f
!***begin prologue     check_norm_2d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for two dimension 
!                      wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_norm_2d_d
  SUBROUTINE check_norm_2d_d(v,norm)
  USE dvrprop_global_it
  USE dvr_shared,           ONLY  : nphy
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1))       :: v
  REAL*8                                   :: ddot
  REAL*8                                   :: norm
!
  norm = ddot(nphy(2)*nphy(1),v,1,v,1)
  write(iout,1) norm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE check_norm_2d_d
!deck check_norm_3d_d.f
!***begin prologue     check_norm_3d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for three dimension wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_norm_3d_d
  SUBROUTINE check_norm_3d_d(v,norm)
  USE dvrprop_global_it
  USE dvr_shared,           ONLY  : nphy
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))        :: v
  REAL*8                                            :: ddot
  REAL*8                                            :: norm
!
  norm = ddot(nphy(3)*nphy(2)*nphy(1),v,1,v,1)
  write(iout,1) norm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE check_norm_3d_d
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
  USE dvrprop_global_it
  USE dvr_global
  USE dvr_shared
  USE plot_wavefunction
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))             :: wave_function
  REAL*8                                 :: norm
  INTEGER                                :: i  
  WRITE(iout,1)
  IF(typ_pak == 'exponential') then
     WRITE(iout,2)
     wave_function = ( grid(1)%pt - x_0(1) )**powr(1)                     &
                                  *                                       &
                      exp( - ( sigma(1) ) * ( grid(1)%pt - x_0(1) ) )
  ELSE IF(typ_pak == 'gaussian') then
     WRITE(iout,3)
     wave_function = ( grid(1)%pt - x_0(1) )**powr(1)                     &
                                  *                                       &
                  exp( - sigma(1) * ( grid(1)%pt - x_0(1) )               &
                                  *                                       &
                                    ( grid(1)%pt - x_0(1) ) )
  END IF
  WRITE(iout,4)
  WRITE(iout,5) sigma(1), x_0(1), powr(1)
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     wave_function  = wave_function * grid(1)%wt
     DO i=1,nphy(1)
        wave_function(i) = wave_function(i) * grid(1)%f(i,i)
     END DO
  END IF
  call check_norm(wave_function,norm)
  wave_function = wave_function * 1.d0/sqrt(norm)  
  title='initial wavepacket'
  CALL print_psi(wave_function)
  write(iout,7) norm
1    FORMAT(/,5X,'initial wavepacket at t=0')
2    FORMAT(/,5X,'the form of the radial packet is:',///,5X,  &
                 'psi = (r - x_0 )**n * exp( - alpha  * ' &
                                           '( r - x_0 ))')
3    FORMAT(/,5X,'the form of the radial packet is:',///,5X,  &
                 'psi = ( r - x_0 )**n * &
                  exp( - alpha * ( r -x_0 ) * ( r - x_0 ) * ' &
                                                            '(r- x_0 ) )')
4    FORMAT(/,1X,'radial wave packet parameters')
5    FORMAT(/,1X,'exponent    = ',e15.8, &
            /,1X,'shift       = ',e15.8, &
            /,1X,'power       = ',i3)
7    FORMAT(/,2x,'Normalization of Initial Wavepacket = ',e15.8,1x,e15.8)
END SUBROUTINE rad_packet
END MODULE normalize

