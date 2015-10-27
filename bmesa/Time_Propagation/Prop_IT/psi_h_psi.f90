!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                               MODULE psi_h_psi
!**begin prologue     psi_h_psi
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       psi_h_psi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                               INTERFACE check_energy
                   MODULE PROCEDURE check_energy_1d_d, check_energy_2d_d, check_energy_3d_d
                          END  INTERFACE check_energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
                               CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck check_energy_1d_d.f
!***begin prologue     check_energy_1d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check energy one dimension wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_energy_1d_d
  SUBROUTINE check_energy_1d_d(v,hv,e)
  USE dvrprop_global_it
  USE dvr_shared,           ONLY  : nphy
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1))             :: v
  REAL*8, DIMENSION(nphy(1))             :: hv
  REAL*8                                 :: ddot
  REAL*8                                 :: e
!
  e = ddot(nphy(1),v,1,hv,1)
  write(iout,1) e
1 format(/,10x,'Energy = ',e20.12)
END SUBROUTINE check_energy_1d_d
!deck check_energy_2d_d.f
!***begin prologue     check_energy_2d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            energy for two dimensional wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_energy_2d_d
  SUBROUTINE check_energy_2d_d(v,hv,g_1,e)
  USE dvrprop_global_it
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1))       :: v
  REAL*8, DIMENSION(nphy(2),nphy(1))       :: hv
  REAL*8, DIMENSION(nphy(1))               :: g_1
  REAL*8                                   :: ddot
  REAL*8                                   :: e, first, last
  INTEGER                                  :: i
!
  e = 0.d0
  DO i=1,nphy(1)
     g_1(i)  =    ddot(nphy(2),v(1,i),1,hv(1,i),1)      
  END DO
  IF( typke == 'fd' ) then
      DO i=1,nphy(1)
         first = v(1,i) * hv(1,i)                       
         last  = v(nphy(2),i) * hv(nphy(2),i)           
         g_1(i) = g_1(i) - .5d0 * ( first + last)
      END DO
      g_1(:) = del * g_1(:)
  END IF
  DO i=1,nphy(1)
     e = e + g_1(i)
  END DO
  IF( typke == 'fd' ) then
      e = e - .5d0 * ( g_1(1) + g_1(nphy(1)) )
      e = del * e
  END IF
  write(iout,1) e
1 format(/,5x,'Energy = ',e20.12)
END SUBROUTINE check_energy_2d_d
!deck check_energy_3d_d.f
!***begin prologue     check_energy_3d_d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            energy for three dimensional wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_energy_3d_d
  SUBROUTINE check_energy_3d_d(v,hv,g_1,g_2,e)
  USE dvrprop_global_it
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))        :: v
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))        :: hv
  REAL*8, DIMENSION(nphy(1))                        :: g_1
  REAL*8, DIMENSION(nphy(2),nphy(1))                :: g_2
  REAL*8                                            :: ddot
  REAL*8                                            :: e, first, last
  INTEGER                                           :: i, j
!
  e = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_2(j,i) =  ddot(nphy(3),v(1,j,i),1,hv(1,j,i),1)                
     END DO
  END DO
  IF( typke == 'fd' ) then
      DO i=1,nphy(1)
         DO j=1,nphy(2)
            first = v(1,j,i) * hv(1,j,i)                                
            last  = v(nphy(3),j,i) * hv(nphy(3),j,i)                    
            g_2(j,i) = g_2(j,i) - .5d0 * ( first + last)
         END DO
      END DO
      g_2(:,:) = del * g_2(:,:)
  END IF
  g_1(:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_1(i) = g_1(i) + g_2(j,i)
     END DO
  END DO
  IF( typke == 'fd' ) then
      g_1(1:nphy(1)) = g_1(1:nphy(1)) - .5d0 * ( g_2(1,1:nphy(1))          &
                                      +                                    &
                           g_2(nphy(2),1:nphy(1)) )
      g_1(:) = del * g_1(:)
  END IF
  DO i =1,nphy(1)
     e = e + g_1(i)
  END DO
  IF( typke == 'fd' ) then
      e = e -.5d0 * ( g_1(1) + g_1(nphy(1)) )
      e = del * e
  END IF
  write(iout,1) e
1 format(/,5x,'Energy = ',e20.12)
END SUBROUTINE check_energy_3d_d
END MODULE psi_h_psi

