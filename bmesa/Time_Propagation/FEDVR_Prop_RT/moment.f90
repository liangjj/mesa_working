!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       MODULE moment
!**begin prologue     moment
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       moment
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      INTERFACE calc_moment
             MODULE PROCEDURE calc_moment_1d_z, calc_moment_2d_z,        &
                              calc_moment_3d_z
                  END INTERFACE calc_moment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck calc_moment_1d_z.f
!***begin prologue     calc_moment_1d_z
!***date written       040706   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate some moments and compare
!                      
!***
!***references
!***routines called
!***end prologue       calc_moment_1d_z
!
  SUBROUTINE calc_moment_1d_z(vec_1d_z,t_calc)
  USE dvrprop_global_rt
  USE dvr_shared
  USE fd_global,           ONLY : del
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(1),2)           :: vec_1d_z
  REAL*8                                 :: t_calc, av, avv, fac_1d_z
  REAL*8                                 :: exav, exavv
  INTEGER                                :: i
!
!
  av = 0.d0
  avv =0.d0
  DO i =1,nphy(1)
     fac_1d_z = vec_1d_z(i,1) * vec_1d_z(i,1)                    &
                          +                                &
              vec_1d_z(i,2) * vec_1d_z(i,2)
     av = av   + grid(1)%pt(i) * fac_1d_z 
     avv = avv + grid(1)%pt(i) * grid(1)%pt(i) * fac_1d_z 
  END DO
  avv = avv - av * av
  avv = SQRT(avv)
  exav = - (beta(1) * t_calc) + x_0(1)
  exavv = .5d0 * ( sigma(1) * sigma(1) +  ( t_calc * t_calc) /  &
                                          ( sigma(1) * sigma(1) ) )
  exavv = SQRT(exavv)
  write(iout,1) t_calc, exav, av, exavv, avv
1 format(/,1x,'Time                      = ',e15.8,        &
         /,1x,'Exact <x>                 = ',e15.8,        &
         /,1x,'Calculated <x>            = ',e15.8,        &
         /,1x,'Exact <x*x> - <x><x>      = ',e15.8         &
         /,1x,'Calculated <x*x> - <x><x> = ',e15.8)
END SUBROUTINE calc_moment_1d_z
!deck calc_moment_2d_z.f
!***begin prologue     calc_moment_2d_z
!***date written       0400706   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate some moments and compare
!                      
!***
!***references
!***routines called
!***end prologue       calc_moment_2d_z
!
  SUBROUTINE calc_moment_2d_z(vec_2d_z,t_calc,g_2)
  USE dvrprop_global_rt
  USE dvr_shared
  USE fd_global,             ONLY  : del
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(2),nphy(1),2)      :: vec_2d_z
  REAL*8, DIMENSION(maxdim,2)               :: g_2
  REAL*8, DIMENSION(2)                      :: av, avv
  REAL*8, DIMENSION(2)                      :: exav, exavv
  REAL*8                                    :: t_calc, f
  REAL*8                                    :: ddot
  INTEGER                                   :: i, j
!
!
  av = 0.d0
  avv = 0.d0
  g_2(1:nphy(1),:) = 0.d0
  DO i =1,nphy(1)
     DO j=1,nphy(2)
        f = ( vec_2d_z(j,i,1) * vec_2d_z(j,i,1)                      &
                            +                                    &
              vec_2d_z(j,i,2) * vec_2d_z(j,i,2) )                    &
                            *                                    &
                     grid(2)%pt(j)
        g_2(i,1) = g_2(i,1) + f 
        g_2(i,2) = g_2(i,2) + f * grid(2)%pt(j)
     END DO
  END DO
  DO i=1,nphy(1)
     av(2) =  av(2)  + g_2(i,1)
     avv(2) = avv(2) + g_2(i,2)
  END DO
  g_2(1:nphy(2),:) = 0.d0
  DO i =1,nphy(2)
     DO j=1,nphy(1)
        f = ( vec_2d_z(i,j,1) * vec_2d_z(i,j,1)                      &
                            +                                    &
              vec_2d_z(i,j,2) * vec_2d_z(i,j,2) )                    &
                            *                                    &
                     grid(1)%pt(j)
        g_2(i,1) = g_2(i,1) + f 
        g_2(i,2) = g_2(i,2) + f * grid(1)%pt(j)
     END DO
  END DO
  DO i=1,nphy(2)
     av(1) =  av(1)  + g_2(i,1)
     avv(1) = avv(1) + g_2(i,2)
  END DO
  avv(1) = avv(1) - av(1) * av(1)
  avv(2) = avv(2) - av(2) * av(2)
  avv(1) = SQRT(avv(1))
  avv(2) = SQRT(avv(2))
  exav(1) = - (beta(1) * t_calc) + x_0(1)
  exav(2) = - (beta(2) * t_calc) + x_0(2)
  exavv(1) = .5d0 * ( sigma(1) * sigma(1) + t_calc * t_calc      &
            /  (sigma(1) * sigma(1)) )
  exavv(2) = .5d0 * ( sigma(2) * sigma(2) + t_calc * t_calc      &
            /  (sigma(2) * sigma(2)) )
  exavv(1) = SQRT(exavv(1))
  exavv(2) = SQRT(exavv(2))
  write(iout,1) t_calc, exav(1), av(1), exavv(1), avv(1)
  write(iout,2) exav(2), av(2), exavv(2), avv(2)
1 format(/,1x,'Time                      = ',e15.8,              &
         /,1x,'Exact <x>                 = ',e15.8,              &
         /,1x,'Calculated <x>            = ',e15.8,              &
         /,1x,'Exact <x*x> - <x><x>      = ',e15.8               &
         /,1x,'Calculated <x*x> - <x><x> = ',e15.8)
2 format(/,1x,'Exact <y>                 = ',e15.8,              &
         /,1x,'Calculated <y>            = ',e15.8,              &
         /,1x,'Exact <y*y> - <y><y>      = ',e15.8               &
         /,1x,'Calculated <y*y> - <y><y> = ',e15.8)
END SUBROUTINE calc_moment_2d_z
!deck calc_moment_3d_z.f
!***begin prologue     calc_moment_3d_z
!***date written       040706   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate some moments and compare
!                      
!***
!***references
!***routines called
!***end prologue       calc_moment_3d_z
!
  SUBROUTINE calc_moment_3d_z(vec_3d_z,t_calc,g_2,g_3)
  USE dvrprop_global_rt
  USE dvr_shared
  USE fd_global,                 ONLY : del
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: vec_3d_z
  REAL*8, DIMENSION(maxdim,2)                    :: g_2
  REAL*8, DIMENSION(maxdim,maxdim,2)             :: g_3
  REAL*8, DIMENSION(3)                           :: av, avv
  REAL*8, DIMENSION(3)                           :: exav, exavv
  REAL*8                                         :: f, ddot
  REAL*8                                         :: t_calc
  INTEGER                                        :: i, j, k
!
!
  av = 0.d0
  avv =0.d0
  g_3(1:nphy(2),1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)   
     DO j =1,nphy(2)
        DO k=1,nphy(3)
           f = ( vec_3d_z(k,j,i,1) * vec_3d_z(k,j,i,1)                     &
                                 +                                     &
                 vec_3d_z(k,j,i,2) * vec_3d_z(k,j,i,2) )                   &
                                 *                                     &
                          grid(3)%pt(k)
           g_3(j,i,1) = g_3(j,i,1) + f
           g_3(j,i,2) = g_3(j,i,2) + f * grid(3)%pt(k)
        END DO
    END DO
  END DO
  g_2(1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_2(i,:) = g_2(i,:) + g_3(j,i,:)
     END DO
  END DO
  DO i=1,nphy(1)
     av(3)   = av(3)  + g_2(i,1)
     avv(3)  = avv(3) + g_2(i,2)
  END DO
  g_3(1:nphy(3),1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)   
     DO j =1,nphy(3)
        DO k=1,nphy(2)
           f = ( vec_3d_z(j,k,i,1) * vec_3d_z(j,k,i,1)                     &
                                 +                                     &
                 vec_3d_z(j,k,i,2) * vec_3d_z(j,k,i,2) )                   &
                                 *                                     &
                          grid(2)%pt(k)
           g_3(j,i,1) = g_3(j,i,1) + f
           g_3(j,i,2) = g_3(j,i,2) + f * grid(2)%pt(k)
        END DO
    END DO
  END DO
  g_2(1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(3)
        g_2(i,:) = g_2(i,:) + g_3(j,i,:)
     END DO
  END DO
  DO i=1,nphy(1)
     av(2)   = av(2)  + g_2(i,1)
     avv(2)  = avv(2) + g_2(i,2)
  END DO
  g_3(1:nphy(3),1:nphy(2),:) = 0.d0
  DO i=1,nphy(3)   
     DO j =1,nphy(2)
        DO k=1,nphy(1)
           f = ( vec_3d_z(i,j,k,1) * vec_3d_z(i,j,k,1)                     &
                                 +                                     &
                 vec_3d_z(i,j,k,2) * vec_3d_z(i,j,k,2) )                   &
                                 *                                     &
                          grid(1)%pt(k)
           g_3(j,i,1) = g_3(j,i,1) + f
           g_3(j,i,2) = g_3(j,i,2) + f * grid(1)%pt(k)
        END DO
    END DO
  END DO
  g_2(1:nphy(3),:) = 0.d0
  DO i=1,nphy(3)
     DO j=1,nphy(2)
        g_2(i,:) = g_2(i,:) + g_3(j,i,:)
     END DO
  END DO
  DO i=1,nphy(1)
     av(1)   = av(1)  + g_2(i,1)
     avv(1)  = avv(1) + g_2(i,2)
  END DO
  avv(1) = avv(1) - av(1) * av(1)
  avv(2) = avv(2) - av(2) * av(2)
  avv(3) = avv(3) - av(3) * av(3)
  avv = abs(avv)
  avv(1) = SQRT(avv(1))
  avv(2) = SQRT(avv(2))
  avv(3) = SQRT(avv(3))
  exav(1) = - (beta(1) * t_calc) + x_0(1)
  exav(2) = - (beta(2) * t_calc) + x_0(2)
  exav(3) = - (beta(3) * t_calc) + x_0(3)
  exavv(1) = .5d0 * ( sigma(1) * sigma(1) + t_calc * t_calc            &
               /  (sigma(1) * sigma(1)) )
  exavv(2) = .5d0 * ( sigma(2) * sigma(2) + t_calc * t_calc            &
               /  (sigma(2) * sigma(2)) )
  exavv(3) = .5d0 * ( sigma(3) * sigma(3) + t_calc * t_calc            &
               /  ( sigma(3) * sigma(3)) )
  exavv(1) = SQRT(exavv(1))
  exavv(2) = SQRT(exavv(2))
  exavv(3) = SQRT(exavv(3))
  write(iout,1) t_calc, exav(1), av(1), exavv(1), avv(1)
  write(iout,2) exav(2), av(2), exavv(2), avv(2)
  write(iout,3) exav(3), av(3), exavv(3), avv(3)
1 format(/,1x,'Time                      = ',e15.8,                    &
         /,1x,'Exact <x>                 = ',e15.8,                    &
         /,1x,'Calculated <x>            = ',e15.8,                    &
         /,1x,'Exact <x*x> - <x><x>      = ',e15.8                     &
         /,1x,'Calculated <x*x> - <x><x> = ',e15.8)
2 format(/,1x,'Exact <y>                 = ',e15.8,                    &
         /,1x,'Calculated <y>            = ',e15.8,                    &
         /,1x,'Exact <y*y> - <y><y>      = ',e15.8                     &
         /,1x,'Calculated <y*y> - <y><y> = ',e15.8)
3 format(/,1x,'Exact <z>                 = ',e15.8,                    &
         /,1x,'Calculated <z>            = ',e15.8,                    &
         /,1x,'Exact <z*z> - <z><z>      = ',e15.8                     &
         /,1x,'Calculated <z*z> - <z><z> = ',e15.8)
END SUBROUTINE calc_moment_3d_z
!deck moment_data
!***begin prologue     moment_data
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            input data for moment calculation.
!                      
!***
!***references
!***routines called
!***end prologue       moment_data
!
  SUBROUTINE moment_data
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE                          
  INTEGER                                :: intkey
  REAL*8                                 :: fpkey
  LOGICAL                                :: dollar
  CHARACTER(LEN=80)                      :: chrkey
!
!
  IF ( dollar('$initial-state',card,title,inp) ) THEN
       CALL fparr(card,'alpha',alpha,spdim,' ')
       CALL fparr(card,'sigma',sigma,spdim,' ')
       CALL fparr(card,'x_0',x_0,spdim,' ')
       CALL fparr(card,'beta',beta,spdim,' ')
       powr=intkey(card,'power',1,' ')
       typ_pak=chrkey(card,'type-radial-packet', &
                           'exponential',' ')
  ELSE
       write(iout,1)
       call lnkerr
  END IF
  sigma = sigma/sqrt(alpha)
  alpha = 1.d0
1 format(/,1x,'This is a test for a separable free wavepacket', &
         /,1x,'Current calculation does not qualify')
END SUBROUTINE moment_data
END MODULE moment

