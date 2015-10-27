!deck mom_2d.f
!***begin prologue     mom_2d
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
!***end prologue       mom_2d
!
  SUBROUTINE mom_2d(vec_2d,fac_2d,av,avv)
  USE dvrprop_global
  USE dvr_shared
!  USE dvr_global
  USE fd_global,             ONLY  : del
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(2),nphy(1),2)      :: vec_2d
  REAL*8, DIMENSION(maxdim,2)               :: fac_2d
  REAL*8, DIMENSION(2)                      :: av, avv
  REAL*8                                    :: f
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
        f = ( vec_2d(j,i,1) * vec_2d(j,i,1)    &
                            +                  &
              vec_2d(j,i,2) * vec_2d(j,i,2) )  &
                            *                  &
                     grid(2)%pt(j)
        fac_2d(i,1) = fac_2d(i,1) + f 
        fac_2d(i,2) = fac_2d(i,2) + f * grid(2)%pt(j)
     END DO
  END DO
  IF (typke == 'fd') THEN
      DO i=1,nphy(1)
         first = ( vec_2d(1,i,1) * vec_2d(1,i,1)    &
                                 +                  &
                   vec_2d(1,i,2) * vec_2d(1,i,2) )  &
                                 *                  &  
                          grid(2)%pt(1)       
         last =  ( vec_2d(nphy(2),i,1) * vec_2d(nphy(2),i,1)      &
                                       +                          &
                   vec_2d(nphy(2),i,2) * vec_2d(nphy(2),i,2) )    &
                                       *                          &
                                grid(2)%pt(nphy(2))
         fac_2d(i,1) =  fac_2d(i,1) - .5d0 * ( first + last )
         first = first * grid(2)%pt(1)
         last  = last  * grid(2)%pt(nphy(2))
         fac_2d(i,2) = fac_2d(i,2)- .5d0 * ( first + last )
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
        f = ( vec_2d(i,j,1) * vec_2d(i,j,1)    &
                            +                  &
              vec_2d(i,j,2) * vec_2d(i,j,2) )  &
                            *                  &
                     grid(1)%pt(j)
        fac_2d(i,1) = fac_2d(i,1) + f 
        fac_2d(i,2) = fac_2d(i,2) + f * grid(1)%pt(j)
     END DO
  END DO
  IF (typke == 'fd') THEN
      DO i=1,nphy(2)
         first = ( vec_2d(i,1,1) * vec_2d(i,1,1)    &
                                 +                  &
                   vec_2d(i,1,2) * vec_2d(i,1,2) )  &
                                 *                  &  
                          grid(1)%pt(1)       
         last =  ( vec_2d(i,nphy(1),1) * vec_2d(i,nphy(1),1)      &
                                       +                          &
                   vec_2d(i,nphy(1),2) * vec_2d(i,nphy(1),2) )    &
                                       *                          &
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
END SUBROUTINE mom_2d
