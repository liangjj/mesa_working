!deck mom_1d.f
!***begin prologue     mom_1d
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
!***end prologue       mom_1d
!
  SUBROUTINE mom_1d(vec_1d,av,avv)
  USE dvrprop_global
  USE dvr_shared
!  USE dvr_global
  USE fd_global,           ONLY : del
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(1),2)           :: vec_1d
  REAL*8                                 :: av, avv, fac_1d
  REAL*8                                 :: first, last
  INTEGER                                :: i
!
!
  av = 0.d0
  avv =0.d0
  DO i =1,nphy(1)
     fac_1d = vec_1d(i,1) * vec_1d(i,1) + vec_1d(i,2) * vec_1d(i,2)
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
END SUBROUTINE mom_1d
