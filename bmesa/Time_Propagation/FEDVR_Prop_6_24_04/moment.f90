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
  IMPLICIT NONE
  CONTAINS
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
  USE dvrprop_global
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
!deck moment_1d.f
!***begin prologue     moment_1d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate some moments and compare
!                      
!***
!***references
!***routines called
!***end prologue       moment_1d
!
  SUBROUTINE moment_1d(vec_1d,t_calc)
  USE dvrprop_global
  USE dvr_shared
  USE fd_global,           ONLY : del
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(1),2)           :: vec_1d
  REAL*8                                 :: t_calc, av, avv, fac_1d
  REAL*8                                 :: exav, exavv
  REAL*8                                 :: first, last
  INTEGER                                :: i
!
!
  av = 0.d0
  avv =0.d0
  DO i =1,nphy(1)
     fac_1d = vec_1d(i,1) * vec_1d(i,1)                    &
                          +                                &
              vec_1d(i,2) * vec_1d(i,2)
     av = av   + grid(1)%pt(i) * fac_1d 
     avv = avv + grid(1)%pt(i) * grid(1)%pt(i) * fac_1d 
  END DO
  IF (typke == 'fd') THEN
      first = ( vec_1d(1,1) * vec_1d(1,1)                  &
                            +                              &
               vec_1d(1,2)  * vec_1d(1,2) )                &
                            *                              & 
                     grid(1)%pt(1)
      last  = ( vec_1d(nphy(1),1) * vec_1d(nphy(1),1)      &
                                  +                        &
                vec_1d(nphy(1),2) * vec_1d(nphy(1),2) )    &
                                  *                        &
                           grid(1)%pt(nphy(1))
      av = av -.5d0 * ( first + last )
      av = av * del
      first = first * grid(1)%pt(1)
      last = last * grid(1)%pt(nphy(1))
      avv = avv -.5d0 * ( first + last )
      avv = avv * del
  END IF
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
END SUBROUTINE moment_1d
!deck moment_2d.f
!***begin prologue     moment_2d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate some moments and compare
!                      
!***
!***references
!***routines called
!***end prologue       moment_2d
!
  SUBROUTINE moment_2d(vec_2d,fac_2d,t_calc)
  USE dvrprop_global
  USE dvr_shared
  USE fd_global,             ONLY  : del
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(2),nphy(1),2)      :: vec_2d
  REAL*8, DIMENSION(maxdim,2)               :: fac_2d
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
  fac_2d(1:nphy(1),:) = 0.d0
  DO i =1,nphy(1)
     DO j=1,nphy(2)
        f = ( vec_2d(j,i,1) * vec_2d(j,i,1)                      &
                            +                                    &
              vec_2d(j,i,2) * vec_2d(j,i,2) )                    &
                            *                                    &
                     grid(2)%pt(j)
        fac_2d(i,1) = fac_2d(i,1) + f 
        fac_2d(i,2) = fac_2d(i,2) + f * grid(2)%pt(j)
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
         fac_2d(i,1) =  fac_2d(i,1)                              &
                                    -                            &
                        .5d0 * ( first + last )
         first = first * grid(2)%pt(1)
         last  = last  * grid(2)%pt(nphy(2))
         fac_2d(i,2) = fac_2d(i,2)                               &
                                    -                            &
                       .5d0 * ( first + last )               
     END DO
     fac_2d(1:nphy(1),:) = del * fac_2d(1:nphy(1),:)
  END IF
  DO i=1,nphy(1)
     av(2) =  av(2)  + fac_2d(i,1)
     avv(2) = avv(2) + fac_2d(i,2)
  END DO
  IF (typke == 'fd') THEN
      av(2)  = av(2)  - .5d0 * ( fac_2d(1,1) + fac_2d(nphy(1),1) )
      avv(2) = avv(2) - .5d0 * ( fac_2d(1,2) + fac_2d(nphy(1),2) )
      av(2)  = del * av(2)
      avv(2) = del * avv(2)
  END IF
  fac_2d(1:nphy(2),:) = 0.d0
  DO i =1,nphy(2)
     DO j=1,nphy(1)
        f = ( vec_2d(i,j,1) * vec_2d(i,j,1)                      &
                            +                                    &
              vec_2d(i,j,2) * vec_2d(i,j,2) )                    &
                            *                                    &
                     grid(1)%pt(j)
        fac_2d(i,1) = fac_2d(i,1) + f 
        fac_2d(i,2) = fac_2d(i,2) + f * grid(1)%pt(j)
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
         fac_2d(i,1) =  fac_2d(i,1) - .5d0 * ( first + last )
         first = first * grid(1)%pt(1)
         last  = last  * grid(1)%pt(nphy(1))
         fac_2d(i,2) = fac_2d(i,2)- .5d0 * ( first + last )
      END DO  
      fac_2d(1:nphy(2),:) = del * fac_2d(1:nphy(2),:)
  END IF
  DO i=1,nphy(2)
     av(1) =  av(1)  + fac_2d(i,1)
     avv(1) = avv(1) + fac_2d(i,2)
  END DO
  IF (typke == 'fd') THEN
      av(1)  = av(1)  - .5d0 * ( fac_2d(1,1) + fac_2d(nphy(2),1) )
      avv(1) = avv(1) - .5d0 * ( fac_2d(1,2) + fac_2d(nphy(2),2) )
      av(1)  = del * av(1)
      avv(1) = del * avv(1)
  END IF
  avv(1) = avv(1) - av(1) * av(1)
  avv(2) = avv(2) - av(2) * av(2)
  avv(1) = avv(1) - av(1) * av(1)
  avv(2) = avv(2) - av(2) * av(2)
  avv(1) = SQRT(avv(1))
  avv(2) = SQRT(avv(2))
  exav(1) = - (beta(1) * t_calc) + x_0(1)
  exav(2) = - (beta(2) * t_calc) + x_0(2)
  exavv(1) = .5d0 * ( sigma(1) * sigma(1) + t_calc * t_calc              &
            /  (sigma(1) * sigma(1)) )
  exavv(2) = .5d0 * ( sigma(2) * sigma(2) + t_calc * t_calc              &
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
END SUBROUTINE moment_2d
!deck moment_3d.f
!***begin prologue     moment_3d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            calculate some moments and compare
!                      
!***
!***references
!***routines called
!***end prologue       moment_3d
!
  SUBROUTINE moment_3d(vec_3d,fac_3d,v_3d,t_calc)
  USE dvrprop_global
  USE dvr_shared
  USE fd_global,                 ONLY : del
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: vec_3d
  REAL*8, DIMENSION(3)                           :: av, avv
  REAL*8, DIMENSION(3)                           :: exav, exavv
  REAL*8, DIMENSION(maxdim,maxdim,2)             :: fac_3d
  REAL*8, DIMENSION(maxdim,2)                    :: v_3d
  REAL*8                                         :: f, ddot, first, last
  REAL*8                                         :: t_calc
  INTEGER                                        :: i, j, k
!
!
  av = 0.d0
  avv =0.d0
  fac_3d(1:nphy(2),1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)   
     DO j =1,nphy(2)
        DO k=1,nphy(3)
           f = ( vec_3d(k,j,i,1) * vec_3d(k,j,i,1)                  &
                                 +                                  &
                 vec_3d(k,j,i,2) * vec_3d(k,j,i,2) )                &
                                 *                                  &
                          grid(3)%pt(k)
           fac_3d(j,i,1) = fac_3d(j,i,1) + f
           fac_3d(j,i,2) = fac_3d(j,i,2) + f * grid(3)%pt(k)
        END DO
    END DO
  END DO
  IF (typke == 'fd') THEN
     DO i=1,nphy(1)
        DO j=1,nphy(2)
           first = ( vec_3d(1,j,i,1) * vec_3d(1,j,i,1)              &
                                     +                              &
                     vec_3d(1,j,i,2) * vec_3d(1,j,i,2) )            &
                                     *                              & 
                              grid(3)%pt(1)       
           last = ( vec_3d(nphy(3),j,i,1) * vec_3d(nphy(3),j,i,1)   &
                                          +                         &
                    vec_3d(nphy(3),j,i,2) * vec_3d(nphy(3),j,i,2) ) &
                                          *                         &
                                   grid(3)%pt(nphy(3))
           fac_3d(j,i,1) =  fac_3d(j,i,1) - .5d0 * ( first + last )
           first = first * grid(3)%pt(1)
           last  = last  * grid(3)%pt(nphy(3))
           fac_3d(j,i,2) = fac_3d(j,i,2)- .5d0 * ( first + last )
        END DO
     END DO
     fac_3d(1:nphy(2),1:nphy(1),:) = del                            &
                                         *                          &
                                     fac_3d(1:nphy(2),1:nphy(1),:)
  END IF
  v_3d(1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        v_3d(i,:) = v_3d(i,:) + fac_3d(j,i,:)
     END DO
  END DO
  IF (typke == 'fd') THEN 
      DO i=1,nphy(1)
         v_3d(i,:) = v_3d(i,:) - .5d0 * ( fac_3d(nphy(2),i,:)       & 
                               +                                    &
                           fac_3d(1,i,:) )
      END DO
      v_3d(1:nphy(1),:) = del * v_3d(1:nphy(1),:)
  END IF
  DO i=1,nphy(1)
     av(3)   = av(3)  + v_3d(i,1)
     avv(3)  = avv(3) + v_3d(i,2)
  END DO
  IF (typke == 'fd') THEN       
      av(3)    = av(3)  -.5d0 * ( v_3d(1,1) + v_3d(nphy(1),1) )
      avv(3)   = avv(3)  -.5d0 * ( v_3d(1,2) + v_3d(nphy(1),2) )
      av(3)  = del * av(3)
      avv(3) = del * avv(3)
  END IF
  fac_3d(1:nphy(3),1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)   
     DO j =1,nphy(3)
        DO k=1,nphy(2)
           f = ( vec_3d(j,k,i,1) * vec_3d(j,k,i,1)                  &
                                 +                                  &
                 vec_3d(j,k,i,2) * vec_3d(j,k,i,2) )                &
                                 *                                  &
                          grid(2)%pt(k)
           fac_3d(j,i,1) = fac_3d(j,i,1) + f
           fac_3d(j,i,2) = fac_3d(j,i,2) + f * grid(2)%pt(k)
        END DO
    END DO
  END DO
  IF (typke == 'fd') THEN
     DO i=1,nphy(1)
        DO j=1,nphy(3)
           first = ( vec_3d(j,1,i,1) * vec_3d(j,1,i,1)              &
                                     +                              &
                     vec_3d(j,1,i,2) * vec_3d(j,1,i,2) )            &
                                     *                              & 
                              grid(2)%pt(1)       
           last = ( vec_3d(j,nphy(2),i,1) * vec_3d(j,nphy(2),i,1)   &
                                          +                         &
                    vec_3d(j,nphy(2),i,2) * vec_3d(j,nphy(2),i,2) ) &
                                          *                         &
                                   grid(2)%pt(nphy(3))
           fac_3d(j,i,1) =  fac_3d(j,i,1) - .5d0 * ( first + last )
           first = first * grid(2)%pt(1)
           last  = last  * grid(2)%pt(nphy(3))
           fac_3d(j,i,2) = fac_3d(j,i,2)- .5d0 * ( first + last )
        END DO
     END DO
     fac_3d(1:nphy(3),1:nphy(1),:) = del                            &
                                         *                          &
                                     fac_3d(1:nphy(3),1:nphy(1),:)
  END IF
  v_3d(1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(3)
        v_3d(i,:) = v_3d(i,:) + fac_3d(j,i,:)
     END DO
  END DO
  IF (typke == 'fd') THEN 
      DO i=1,nphy(1)
         v_3d(i,:) = v_3d(i,:) - .5d0 * ( fac_3d(nphy(3),i,:)       &
                               +                                    &
                           fac_3d(1,i,:) )
      END DO
      v_3d(1:nphy(1),:) = del * v_3d(1:nphy(1),:)
  END IF
  DO i=1,nphy(1)
     av(2)   = av(2)  + v_3d(i,1)
     avv(2)  = avv(2) + v_3d(i,2)
  END DO
  IF (typke == 'fd') THEN       
      av(2)    = av(2)  -.5d0 * ( v_3d(1,1) + v_3d(nphy(1),1) )
      avv(2)   = avv(2)  -.5d0 * ( v_3d(1,2) + v_3d(nphy(1),2) )
      av(2)  = del * av(2)
      avv(2) = del * avv(2)
  END IF
  fac_3d(1:nphy(3),1:nphy(2),:) = 0.d0
  DO i=1,nphy(3)   
     DO j =1,nphy(2)
        DO k=1,nphy(1)
           f = ( vec_3d(i,j,k,1) * vec_3d(i,j,k,1)                  &
                                 +                                  &
                 vec_3d(i,j,k,2) * vec_3d(i,j,k,2) )                &
                                 *                                  &
                          grid(1)%pt(k)
           fac_3d(j,i,1) = fac_3d(j,i,1) + f
           fac_3d(j,i,2) = fac_3d(j,i,2) + f * grid(1)%pt(k)
        END DO
    END DO
  END DO
  IF (typke == 'fd') THEN
     DO i=1,nphy(3)
        DO j=1,nphy(2)
           first = ( vec_3d(i,j,1,1) * vec_3d(i,j,1,1)              &
                                     +                              &
                     vec_3d(i,j,1,2) * vec_3d(i,j,1,2) )            &
                                     *                              & 
                              grid(1)%pt(1)       
           last = ( vec_3d(i,j,nphy(1),1) * vec_3d(i,j,nphy(1),1)   &
                                          +                         &
                    vec_3d(i,j,nphy(1),2) * vec_3d(i,j,nphy(1),2) ) &
                                          *                         &
                                   grid(2)%pt(nphy(3))
           fac_3d(j,i,1) =  fac_3d(j,i,1) - .5d0 * ( first + last )
           first = first * grid(1)%pt(1)
           last  = last  * grid(1)%pt(nphy(3))
           fac_3d(j,i,2) = fac_3d(j,i,2)- .5d0 * ( first + last )
        END DO
     END DO
     fac_3d(1:nphy(2),1:nphy(3),:) = del                            &
                                         *                          &
                                     fac_3d(1:nphy(2),1:nphy(3),:)
  END IF
  v_3d(1:nphy(3),:) = 0.d0
  DO i=1,nphy(3)
     DO j=1,nphy(2)
        v_3d(i,:) = v_3d(i,:) + fac_3d(j,i,:)
     END DO
  END DO
  IF (typke == 'fd') THEN 
      DO i=1,nphy(1)
         v_3d(i,:) = v_3d(i,:) - .5d0 * ( fac_3d(nphy(2),i,:)       &
                                               +                    &
                                          fac_3d(1,i,:) )
      END DO
      v_3d(1:nphy(3),:) = del * v_3d(1:nphy(3),:)
  END IF
  DO i=1,nphy(1)
     av(1)   = av(1)  + v_3d(i,1)
     avv(1)  = avv(1) + v_3d(i,2)
  END DO
  IF (typke == 'fd') THEN       
      av(1)    = av(1)  -.5d0 * ( v_3d(1,1) + v_3d(nphy(1),1) )
      avv(1)   = avv(1)  -.5d0 * ( v_3d(1,2) + v_3d(nphy(1),2) )
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
  exavv(1) = .5d0 * ( sigma(1) * sigma(1) + t_calc * t_calc              &
               /  (sigma(1) * sigma(1)) )
  exavv(2) = .5d0 * ( sigma(2) * sigma(2) + t_calc * t_calc              &
               /  (sigma(2) * sigma(2)) )
  exavv(3) = .5d0 * ( sigma(3) * sigma(3) + t_calc * t_calc              &
               /  ( sigma(3) * sigma(3)) )
  exavv(1) = SQRT(exavv(1))
  exavv(2) = SQRT(exavv(2))
  exavv(3) = SQRT(exavv(3))
  write(iout,1) t_calc, exav(1), av(1), exavv(1), avv(1)
  write(iout,2) exav(2), av(2), exavv(2), avv(2)
  write(iout,3) exav(3), av(3), exavv(3), avv(3)
1 format(/,1x,'Time                      = ',e15.8,              &
         /,1x,'Exact <x>                 = ',e15.8,              &
         /,1x,'Calculated <x>            = ',e15.8,              &
         /,1x,'Exact <x*x> - <x><x>      = ',e15.8               &
         /,1x,'Calculated <x*x> - <x><x> = ',e15.8)
2 format(/,1x,'Exact <y>                 = ',e15.8,              &
         /,1x,'Calculated <y>            = ',e15.8,              &
         /,1x,'Exact <y*y> - <y><y>      = ',e15.8               &
         /,1x,'Calculated <y*y> - <y><y> = ',e15.8)
3 format(/,1x,'Exact <z>                 = ',e15.8,              &
         /,1x,'Calculated <z>            = ',e15.8,              &
         /,1x,'Exact <z*z> - <z><z>      = ',e15.8               &
         /,1x,'Calculated <z*z> - <z><z> = ',e15.8)
END SUBROUTINE moment_3d
END MODULE moment

