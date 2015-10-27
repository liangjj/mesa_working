!deck check_norm_2d.f
!***begin prologue     check_norm_2d
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
!***end prologue       check_norm_2d
  SUBROUTINE check_norm_2d(v,norm,g_1)
  USE dvrprop_global_rt
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),2)   :: v
  REAL*8, DIMENSION(nphy(1))             :: g_1
  REAL*8                                 :: ddot
  REAL*8                                 :: norm, first, last
  INTEGER                                :: i
!
  norm = 0.d0
  DO i=1,nphy(1)
     g_1(i) =    ddot(nphy(2),v(1,i,1),1,v(1,i,1),1)       &
                                 +                         &
                 ddot(nphy(2),v(1,i,2),1,v(1,i,2),1)
  END DO

  IF( typke == 'fd' ) then
      DO i=1,nphy(1)
         first = v(1,i,1) * v(1,i,1)                       &
                          +                                &
                 v(1,i,2) * v(1,i,2)
         last  = v(nphy(2),i,1) * v(nphy(2),i,1)           &
                                +                          &
                 v(nphy(2),i,2) * v(nphy(2),i,2)
         g_1(i) = g_1(i) - .5d0 * ( first + last)
      END DO
      g_1(:) = del * g_1(:)
      write(iout,*) g_1
  END IF
  DO i=1,nphy(1)
     norm = norm + g_1(i)
  END DO
  IF( typke == 'fd' ) then
      norm = norm - .5d0 * ( g_1(1) + g_1(nphy(1)) )
      norm = del * norm
  END IF
  write(iout,1) norm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE check_norm_2d
