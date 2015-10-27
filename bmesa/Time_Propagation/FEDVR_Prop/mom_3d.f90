!deck mom_3d.f
!***begin prologue     mom_3d
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
!***end prologue       mom_3d
!
  SUBROUTINE mom_3d(vec_3d,fac_3d,v_3d,av,avv)
  USE dvrprop_global
  USE dvr_shared
  USE fd_global,                 ONLY : del
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),2)   :: vec_3d
  REAL*8, DIMENSION(3)                           :: av, avv
  REAL*8, DIMENSION(maxdim,maxdim,2)             :: fac_3d
  REAL*8, DIMENSION(maxdim,2)                    :: v_3d
  REAL*8                                         :: f, ddot, first, last
  INTEGER                                        :: i, j, k
!
!
  av = 0.d0
  avv =0.d0
  fac_3d(1:nphy(2),1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)   
     DO j =1,nphy(2)
        DO k=1,nphy(3)
           f = ( vec_3d(k,j,i,1) * vec_3d(k,j,i,1)    &
                                 +                    &
                 vec_3d(k,j,i,2) * vec_3d(k,j,i,2) )  &
                                 *                    &
                          grid(3)%pt(k)
           fac_3d(j,i,1) = fac_3d(j,i,1) + f
           fac_3d(j,i,2) = fac_3d(j,i,2) + f * grid(3)%pt(k)
        END DO
    END DO
  END DO
  IF (typke == 'fd') THEN
     DO i=1,nphy(1)
        DO j=1,nphy(2)
           first = ( vec_3d(1,j,i,1) * vec_3d(1,j,i,1)    &
                                     +                    &
                     vec_3d(1,j,i,2) * vec_3d(1,j,i,2) )  &
                                     *                    &  
                              grid(3)%pt(1)       
           last = ( vec_3d(nphy(3),j,i,1) * vec_3d(nphy(3),j,i,1)      &
                                          +                            &
                    vec_3d(nphy(3),j,i,2) * vec_3d(nphy(3),j,i,2) )    &
                                          *                            &
                                   grid(3)%pt(nphy(3))
           fac_3d(j,i,1) =  fac_3d(j,i,1) - .5d0 * ( first + last )
           first = first * grid(3)%pt(1)
           last  = last  * grid(3)%pt(nphy(3))
           fac_3d(j,i,2) = fac_3d(j,i,2)- .5d0 * ( first + last )
        END DO
     END DO
     fac_3d(1:nphy(2),1:nphy(1),:) = del *  fac_3d(1:nphy(2),1:nphy(1),:)
  END IF
  v_3d(1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(2)
        v_3d(i,:) = v_3d(i,:) + fac_3d(j,i,:)
     END DO
  END DO
  IF (typke == 'fd') THEN 
      DO i=1,nphy(1)
         v_3d(i,:) = v_3d(i,:) - .5d0 * ( fac_3d(nphy(2),i,:) + fac_3d(1,i,:) )
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
           f = ( vec_3d(j,k,i,1) * vec_3d(j,k,i,1)    &
                                 +                    &
                 vec_3d(j,k,i,2) * vec_3d(j,k,i,2) )  &
                                 *                    &
                          grid(2)%pt(k)
           fac_3d(j,i,1) = fac_3d(j,i,1) + f
           fac_3d(j,i,2) = fac_3d(j,i,2) + f * grid(2)%pt(k)
        END DO
    END DO
  END DO
  IF (typke == 'fd') THEN
     DO i=1,nphy(1)
        DO j=1,nphy(3)
           first = ( vec_3d(j,1,i,1) * vec_3d(j,1,i,1)    &
                                     +                    &
                     vec_3d(j,1,i,2) * vec_3d(j,1,i,2) )  &
                                     *                    &  
                              grid(2)%pt(1)       
           last = ( vec_3d(j,nphy(2),i,1) * vec_3d(j,nphy(2),i,1)      &
                                          +                            &
                    vec_3d(j,nphy(2),i,2) * vec_3d(j,nphy(2),i,2) )    &
                                          *                            &
                                   grid(2)%pt(nphy(3))
           fac_3d(j,i,1) =  fac_3d(j,i,1) - .5d0 * ( first + last )
           first = first * grid(2)%pt(1)
           last  = last  * grid(2)%pt(nphy(3))
           fac_3d(j,i,2) = fac_3d(j,i,2)- .5d0 * ( first + last )
        END DO
     END DO
     fac_3d(1:nphy(3),1:nphy(1),:) = del *  fac_3d(1:nphy(3),1:nphy(1),:)
  END IF
  v_3d(1:nphy(1),:) = 0.d0
  DO i=1,nphy(1)
     DO j=1,nphy(3)
        v_3d(i,:) = v_3d(i,:) + fac_3d(j,i,:)
     END DO
  END DO
  IF (typke == 'fd') THEN 
      DO i=1,nphy(1)
         v_3d(i,:) = v_3d(i,:) - .5d0 * ( fac_3d(nphy(3),i,:) + fac_3d(1,i,:) )
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
           f = ( vec_3d(i,j,k,1) * vec_3d(i,j,k,1)    &
                                 +                    &
                 vec_3d(i,j,k,2) * vec_3d(i,j,k,2) )  &
                                 *                    &
                          grid(1)%pt(k)
           fac_3d(j,i,1) = fac_3d(j,i,1) + f
           fac_3d(j,i,2) = fac_3d(j,i,2) + f * grid(1)%pt(k)
        END DO
    END DO
  END DO
  IF (typke == 'fd') THEN
     DO i=1,nphy(3)
        DO j=1,nphy(2)
           first = ( vec_3d(i,j,1,1) * vec_3d(i,j,1,1)    &
                                     +                    &
                     vec_3d(i,j,1,2) * vec_3d(i,j,1,2) )  &
                                     *                    &  
                              grid(1)%pt(1)       
           last = ( vec_3d(i,j,nphy(1),1) * vec_3d(i,j,nphy(1),1)      &
                                          +                            &
                    vec_3d(i,j,nphy(1),2) * vec_3d(i,j,nphy(1),2) )    &
                                          *                            &
                                   grid(2)%pt(nphy(3))
           fac_3d(j,i,1) =  fac_3d(j,i,1) - .5d0 * ( first + last )
           first = first * grid(1)%pt(1)
           last  = last  * grid(1)%pt(nphy(3))
           fac_3d(j,i,2) = fac_3d(j,i,2)- .5d0 * ( first + last )
        END DO
     END DO
     fac_3d(1:nphy(2),1:nphy(3),:) = del *  fac_3d(1:nphy(2),1:nphy(3),:)
  END IF
  v_3d(1:nphy(3),:) = 0.d0
  DO i=1,nphy(3)
     DO j=1,nphy(2)
        v_3d(i,:) = v_3d(i,:) + fac_3d(j,i,:)
     END DO
  END DO
  IF (typke == 'fd') THEN 
      DO i=1,nphy(1)
         v_3d(i,:) = v_3d(i,:) - .5d0 * ( fac_3d(nphy(2),i,:) + fac_3d(1,i,:) )
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
END SUBROUTINE mom_3d
