!deck calc_moment_2d.f
!***begin prologue     calc_moment_2d
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
!***end prologue       calc_moment_2d
!
  SUBROUTINE calc_moment_2d(vec_2d,t_calc,g_2)
  USE dvrprop_global_rt
  USE dvr_shared
  USE fd_global,             ONLY  : del
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(2),nphy(1),2)      :: vec_2d
  REAL*8, DIMENSION(maxdim,2)               :: g_2
  REAL*8, DIMENSION(2)                      :: av, avv
  REAL*8, DIMENSION(2)                      :: exav, exavv
  REAL*8                                    :: t_calc, f
  REAL*8                                    :: ddot
  REAL*8                                    :: first, last
  INTEGER                                   :: i, j
!
!
  av = 0.d0
  avv = 0.d0
  g_2(1:nphy(1),:) = 0.d0
  DO i =1,nphy(1)
     DO j=1,nphy(2)
        f = ( vec_2d(j,i,1) * vec_2d(j,i,1)                      &
                            +                                    &
              vec_2d(j,i,2) * vec_2d(j,i,2) )                    &
                            *                                    &
                     grid(2)%pt(j)
        g_2(i,1) = g_2(i,1) + f 
        g_2(i,2) = g_2(i,2) + f * grid(2)%pt(j)
     END DO
  END DO
  IF (typke == 'fd') THEN
      DO i=1,nphy(1)
         first = ( vec_2d(1,i,1) * vec_2d(1,i,1)                 &
                                 +                               &
                   vec_2d(1,i,2) * vec_2d(1,i,2) )               &
                                 *                               &  
                          grid(2)%pt(1)       
         last =  ( vec_2d(nphy(2),i,1) * vec_2d(nphy(2),i,1)     &
                                       +                         &
                   vec_2d(nphy(2),i,2) * vec_2d(nphy(2),i,2) )   &
                                       *                         &
                                grid(2)%pt(nphy(2))
         g_2(i,1) =  g_2(i,1)                                    &
                                    -                            &
                        .5d0 * ( first + last )
         first = first * grid(2)%pt(1)
         last  = last  * grid(2)%pt(nphy(2))
         g_2(i,2) = g_2(i,2)                                     &
                                    -                            &
                       .5d0 * ( first + last )               
     END DO
     g_2(1:nphy(1),:) = del * g_2(1:nphy(1),:)
  END IF
  DO i=1,nphy(1)
     av(2) =  av(2)  + g_2(i,1)
     avv(2) = avv(2) + g_2(i,2)
  END DO
  IF (typke == 'fd') THEN
      av(2)  = av(2)  - .5d0 * ( g_2(1,1) + g_2(nphy(1),1) )
      avv(2) = avv(2) - .5d0 * ( g_2(1,2) + g_2(nphy(1),2) )
      av(2)  = del * av(2)
      avv(2) = del * avv(2)
  END IF
  g_2(1:nphy(2),:) = 0.d0
  DO i =1,nphy(2)
     DO j=1,nphy(1)
        f = ( vec_2d(i,j,1) * vec_2d(i,j,1)                      &
                            +                                    &
              vec_2d(i,j,2) * vec_2d(i,j,2) )                    &
                            *                                    &
                     grid(1)%pt(j)
        g_2(i,1) = g_2(i,1) + f 
        g_2(i,2) = g_2(i,2) + f * grid(1)%pt(j)
     END DO
  END DO
  IF (typke == 'fd') THEN
      DO i=1,nphy(2)
         first = ( vec_2d(i,1,1) * vec_2d(i,1,1)                 &
                                 +                               &
                   vec_2d(i,1,2) * vec_2d(i,1,2) )               &
                                 *                               &  
                          grid(1)%pt(1)       
         last =  ( vec_2d(i,nphy(1),1) * vec_2d(i,nphy(1),1)     &
                                       +                         &
                   vec_2d(i,nphy(1),2) * vec_2d(i,nphy(1),2) )   &
                                       *                         &
                                grid(1)%pt(nphy(1))
         g_2(i,1) =  g_2(i,1) - .5d0 * ( first + last )
         first = first * grid(1)%pt(1)
         last  = last  * grid(1)%pt(nphy(1))
         g_2(i,2) = g_2(i,2)- .5d0 * ( first + last )
      END DO  
      g_2(1:nphy(2),:) = del * g_2(1:nphy(2),:)
  END IF
  DO i=1,nphy(2)
     av(1) =  av(1)  + g_2(i,1)
     avv(1) = avv(1) + g_2(i,2)
  END DO
  IF (typke == 'fd') THEN
      av(1)  = av(1)  - .5d0 * ( g_2(1,1) + g_2(nphy(2),1) )
      avv(1) = avv(1) - .5d0 * ( g_2(1,2) + g_2(nphy(2),2) )
      av(1)  = del * av(1)
      avv(1) = del * avv(1)
  END IF
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
END SUBROUTINE calc_moment_2d
