!deck calc_moment_1d.f
!***begin prologue     calc_moment_1d
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
!***end prologue       calc_moment_1d
!
  SUBROUTINE calc_moment_1d(vec_1d,t_calc)
  USE dvrprop_global_rt
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
END SUBROUTINE calc_moment_1d
