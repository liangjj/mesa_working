!deck calc_moment_3d.f
!***begin prologue     calc_moment_3d
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
!***end prologue       calc_moment_3d
!
  SUBROUTINE calc_moment_3d(vec_3d,t_calc,g_2,g_3)
  USE dvrprop_global_rt
  USE dvr_shared
  USE fd_global,                 ONLY : del
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: vec_3d
  REAL*8, DIMENSION(maxdim,2)                    :: g_2
  REAL*8, DIMENSION(maxdim,maxdim,2)             :: g_3
  REAL*8, DIMENSION(3)                           :: av, avv
  REAL*8, DIMENSION(3)                           :: exav, exavv
  REAL*8                                         :: f, ddot, first, last
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
           f = ( vec_3d(k,j,i,1) * vec_3d(k,j,i,1)                     &
                                 +                                     &
                 vec_3d(k,j,i,2) * vec_3d(k,j,i,2) )                   &
                                 *                                     &
                          grid(3)%pt(k)
           g_3(j,i,1) = g_3(j,i,1) + f
           g_3(j,i,2) = g_3(j,i,2) + f * grid(3)%pt(k)
        END DO
    END DO
  END DO
  IF (typke == 'fd') THEN
     DO i=1,nphy(1)
        DO j=1,nphy(2)
           first = ( vec_3d(1,j,i,1) * vec_3d(1,j,i,1)                 &
                                     +                                 &
                     vec_3d(1,j,i,2) * vec_3d(1,j,i,2) )               &
                                     *                                 & 
                              grid(3)%pt(1)       
           last = ( vec_3d(nphy(3),j,i,1) * vec_3d(nphy(3),j,i,1)      &
                                          +                            &
                    vec_3d(nphy(3),j,i,2) * vec_3d(nphy(3),j,i,2) )    &
                                          *                            &
                                   grid(3)%pt(nphy(3))
           g_3(j,i,1) =  g_3(j,i,1) - .5d0 * ( first + last )
           first = first * grid(3)%pt(1)
           last  = last  * grid(3)%pt(nphy(3))
           g_3(j,i,2) = g_3(j,i,2)- .5d0 * ( first + last )
        END DO
     END DO
     g_3(1:nphy(2),1:nphy(1),:) = del                                  &
                                         *                             &
                                    g_3(1:nphy(2),1:nphy(1),:)
  END IF
  g_2(1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        g_2(i,:) = g_2(i,:) + g_3(j,i,:)
     END DO
  END DO
  IF (typke == 'fd') THEN 
      DO i=1,nphy(1)
         g_2(i,:) = g_2(i,:) - .5d0 * ( g_3(nphy(2),i,:)               & 
                                       +                               &
                                 g_3(1,i,:) )
      END DO
      g_2(1:nphy(1),:) = del * g_2(1:nphy(1),:)
  END IF
  DO i=1,nphy(1)
     av(3)   = av(3)  + g_2(i,1)
     avv(3)  = avv(3) + g_2(i,2)
  END DO
  IF (typke == 'fd') THEN       
      av(3)    = av(3)  -.5d0 * ( g_2(1,1) + g_2(nphy(1),1) )
      avv(3)   = avv(3)  -.5d0 * ( g_2(1,2) + g_2(nphy(1),2) )
      av(3)  = del * av(3)
      avv(3) = del * avv(3)
  END IF
  g_3(1:nphy(3),1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)   
     DO j =1,nphy(3)
        DO k=1,nphy(2)
           f = ( vec_3d(j,k,i,1) * vec_3d(j,k,i,1)                     &
                                 +                                     &
                 vec_3d(j,k,i,2) * vec_3d(j,k,i,2) )                   &
                                 *                                     &
                          grid(2)%pt(k)
           g_3(j,i,1) = g_3(j,i,1) + f
           g_3(j,i,2) = g_3(j,i,2) + f * grid(2)%pt(k)
        END DO
    END DO
  END DO
  IF (typke == 'fd') THEN
     DO i=1,nphy(1)
        DO j=1,nphy(3)
           first = ( vec_3d(j,1,i,1) * vec_3d(j,1,i,1)                 &
                                     +                                 &
                     vec_3d(j,1,i,2) * vec_3d(j,1,i,2) )               &
                                     *                                 & 
                              grid(2)%pt(1)       
           last = ( vec_3d(j,nphy(2),i,1) * vec_3d(j,nphy(2),i,1)      &
                                          +                            &
                    vec_3d(j,nphy(2),i,2) * vec_3d(j,nphy(2),i,2) )    &
                                          *                            &
                                   grid(2)%pt(nphy(3))
           g_3(j,i,1) =  g_3(j,i,1) - .5d0 * ( first + last )
           first = first * grid(2)%pt(1)
           last  = last  * grid(2)%pt(nphy(3))
           g_3(j,i,2) = g_3(j,i,2)- .5d0 * ( first + last )
        END DO
     END DO
     g_3(1:nphy(3),1:nphy(1),:) = del                                  &
                                         *                             &
                                g_3(1:nphy(3),1:nphy(1),:)
  END IF
  g_2(1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(3)
        g_2(i,:) = g_2(i,:) + g_3(j,i,:)
     END DO
  END DO
  IF (typke == 'fd') THEN 
      DO i=1,nphy(1)
         g_2(i,:) = g_2(i,:) - .5d0 * ( g_3(nphy(3),i,:)               &
                                       +                               &
                                  g_3(1,i,:) )
      END DO
      g_2(1:nphy(1),:) = del * g_2(1:nphy(1),:)
  END IF
  DO i=1,nphy(1)
     av(2)   = av(2)  + g_2(i,1)
     avv(2)  = avv(2) + g_2(i,2)
  END DO
  IF (typke == 'fd') THEN       
      av(2)    = av(2)  -.5d0 * ( g_2(1,1) + g_2(nphy(1),1) )
      avv(2)   = avv(2)  -.5d0 * ( g_2(1,2) + g_2(nphy(1),2) )
      av(2)  = del * av(2)
      avv(2) = del * avv(2)
  END IF
  g_3(1:nphy(3),1:nphy(2),:) = 0.d0
  DO i=1,nphy(3)   
     DO j =1,nphy(2)
        DO k=1,nphy(1)
           f = ( vec_3d(i,j,k,1) * vec_3d(i,j,k,1)                     &
                                 +                                     &
                 vec_3d(i,j,k,2) * vec_3d(i,j,k,2) )                   &
                                 *                                     &
                          grid(1)%pt(k)
           g_3(j,i,1) = g_3(j,i,1) + f
           g_3(j,i,2) = g_3(j,i,2) + f * grid(1)%pt(k)
        END DO
    END DO
  END DO
  IF (typke == 'fd') THEN
     DO i=1,nphy(3)
        DO j=1,nphy(2)
           first = ( vec_3d(i,j,1,1) * vec_3d(i,j,1,1)                 &
                                     +                                 &
                     vec_3d(i,j,1,2) * vec_3d(i,j,1,2) )               &
                                     *                                 & 
                              grid(1)%pt(1)       
           last = ( vec_3d(i,j,nphy(1),1) * vec_3d(i,j,nphy(1),1)      &
                                          +                            &
                    vec_3d(i,j,nphy(1),2) * vec_3d(i,j,nphy(1),2) )    &
                                          *                            &
                                   grid(2)%pt(nphy(3))
           g_3(j,i,1) =  g_3(j,i,1) - .5d0 * ( first + last )
           first = first * grid(1)%pt(1)
           last  = last  * grid(1)%pt(nphy(3))
           g_3(j,i,2) = g_3(j,i,2)- .5d0 * ( first + last )
        END DO
     END DO
     g_3(1:nphy(2),1:nphy(3),:) = del                                  &
                                         *                             &
                                 g_3(1:nphy(2),1:nphy(3),:)
  END IF
  g_2(1:nphy(3),:) = 0.d0
  DO i=1,nphy(3)
     DO j=1,nphy(2)
        g_2(i,:) = g_2(i,:) + g_3(j,i,:)
     END DO
  END DO
  IF (typke == 'fd') THEN 
      DO i=1,nphy(1)
         g_2(i,:) = g_2(i,:) - .5d0 * ( g_3(nphy(2),i,:)               &
                                               +                       &
                                          g_3(1,i,:) )
      END DO
      g_2(1:nphy(3),:) = del * g_2(1:nphy(3),:)
  END IF
  DO i=1,nphy(1)
     av(1)   = av(1)  + g_2(i,1)
     avv(1)  = avv(1) + g_2(i,2)
  END DO
  IF (typke == 'fd') THEN       
      av(1)    = av(1)  -.5d0 * ( g_2(1,1) + g_2(nphy(1),1) )
      avv(1)   = avv(1)  -.5d0 * ( g_2(1,2) + g_2(nphy(1),2) )
      av(1)  = del * av(1)
      avv(1) = del * avv(1)
  END IF
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
END SUBROUTINE calc_moment_3d
