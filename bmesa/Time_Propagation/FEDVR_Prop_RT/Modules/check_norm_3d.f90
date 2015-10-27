!deck check_norm_3d.f
!***begin prologue     check_norm_3d
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            check normalization for three dimension wavefunction
!***
!***references
!***routines called    ddot
!***end prologue       check_norm_3d
  SUBROUTINE check_norm_3d(v,norm,g_1,g_2)
  USE dvrprop_global_rt
  USE fd_global,            ONLY  : del
  USE dvr_shared,           ONLY  : nphy, typke
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: v
  REAL*8, DIMENSION(nphy(1))                     :: g_1
  REAL*8, DIMENSION(nphy(2),nphy(1))             :: g_2
  REAL*8                                         :: ddot
  REAL*8                                         :: norm, first, last
  INTEGER                                        :: i, j
!
  norm = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_2(j,i) =  ddot(nphy(3),v(1,j,i,1),1,v(1,j,i,1),1)                &
                                          +                                   &
                       ddot(nphy(3),v(1,j,i,2),1,v(1,j,i,2),1)
     END DO
  END DO
  IF( typke == 'fd' ) then
      DO i=1,nphy(1)
         DO j=1,nphy(2)
            first = v(1,j,i,1) * v(1,j,i,1)                                   &
                               +                                              &
                    v(1,j,i,2) * v(1,j,i,2)
            last  = v(nphy(3),j,i,1) * v(nphy(3),j,i,1)                       &
                                    +                                         &
                    v(nphy(3),j,i,2) * v(nphy(3),j,i,2)
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
      g_1(1:nphy(1)) = g_1(1:nphy(1)) - .5d0 * ( g_2(1,1:nphy(1))             &
                                            +                                 &
                                 g_2(nphy(2),1:nphy(1)) )
      g_1(:) = del * g_1(:)
  END IF
  DO i =1,nphy(1)
     norm = norm + g_1(i)
  END DO
  IF( typke == 'fd' ) then
      norm = norm -.5d0 * ( g_1(1) + g_1(nphy(1)) )
      norm = del * norm
  END IF
  write(iout,1) norm
1 format(/,5x,'Normalization Integral = ',e20.12)
END SUBROUTINE check_norm_3d
