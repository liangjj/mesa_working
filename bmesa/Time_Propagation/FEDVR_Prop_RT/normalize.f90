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
                        INTERFACE check_norm
             MODULE PROCEDURE check_norm_1d_z,check_norm_2d_z,    &
                              check_norm_3d_z
                    END INTERFACE check_norm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck check_norm_1d.f
!***begin prologue     check_norm_1d_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for one dimension wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_norm_1d_z
  SUBROUTINE check_norm_1d_z(v,norm)
  USE dvrprop_global_rt
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),2)           :: v
  REAL*8                                 :: ddot
  REAL*8                                 :: norm, first, last
!
  norm = ddot(nphy(1),v(1,1),1,v(1,1),1)               &
                         +                             &
         ddot(nphy(1),v(1,2),1,v(1,2),1)
  write(iout,1) norm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE check_norm_1d_z
!deck check_norm_2d_z.f
!***begin prologue     check_norm_2d_z
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
!***end prologue       check_norm_2d_z
  SUBROUTINE check_norm_2d_z(v,norm,g_1)
  USE dvrprop_global_rt
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: v
  REAL*8, DIMENSION(nphy(1))             :: g_1
  REAL*8                                 :: ddot
  REAL*8                                 :: norm
  INTEGER                                :: i
!
  norm = 0.d0
  DO i=1,nphy(1)
     g_1(i) =    ddot(nphy(2),v(1,i,1),1,v(1,i,1),1)       &
                                 +                         &
                 ddot(nphy(2),v(1,i,2),1,v(1,i,2),1)
  END DO
  DO i=1,nphy(1)
     norm = norm + g_1(i)
  END DO
  write(iout,1) norm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE check_norm_2d_z
!deck check_norm_3d_z.f
!***begin prologue     check_norm_3d_z
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for three dimension wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_norm_3d_z
  SUBROUTINE check_norm_3d_z(v,norm,g_1,g_2)
  USE dvrprop_global_rt
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: v
  REAL*8, DIMENSION(nphy(1))                     :: g_1
  REAL*8, DIMENSION(nphy(2),nphy(1))             :: g_2
  REAL*8                                         :: ddot
  REAL*8                                         :: norm
  INTEGER                                        :: i, j
!
  norm = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_2(j,i) =  ddot(nphy(3),v(1,j,i,1),1,v(1,j,i,1),1)                &
                                          +                                &
                       ddot(nphy(3),v(1,j,i,2),1,v(1,j,i,2),1)
     END DO
  END DO
  g_1(:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_1(i) = g_1(i) + g_2(j,i)
     END DO
  END DO
  DO i =1,nphy(1)
     norm = norm + g_1(i)
  END DO
  write(iout,1) norm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE check_norm_3d_z
END MODULE normalize
