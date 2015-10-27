!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                       MODULE moment
                   USE dvrprop_global
                   USE dvr_global
                   USE dvr_shared
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
                      INTERFACE calculate_moment
             MODULE PROCEDURE calculate_moment_d,        &
                              calculate_moment_z
                  END INTERFACE calculate_moment
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck calculate_moment_d.f
!***begin prologue     calculate_moment_d
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
!***end prologue       calculate_moment_d
!
  SUBROUTINE calculate_moment_d(vec_d,t_calc)
  IMPLICIT NONE                          
  REAL*8, DIMENSION(n3d)                 :: vec_d
  REAL*8                                 :: t_calc
  IF (spdim == 1 ) THEN
      CALL calculate_moment_1d_d(vec_d,t_calc)
  ELSE IF (spdim == 2 ) THEN
      CALL calculate_moment_2d_d(vec_d,t_calc)
  ELSE IF ( spdim == 3 ) THEN
      CALL calculate_moment_3d_d(vec_d,t_calc)
  END IF
END SUBROUTINE calculate_moment_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck calculate_moment_z.f
!***begin prologue     calculate_moment_z
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
!***end prologue       calculate_moment_z
!
  SUBROUTINE calculate_moment_z(vec_z,t_calc)
  IMPLICIT NONE                          
  COMPLEX*16, DIMENSION(n3d)             :: vec_z
  REAL*8                                 :: t_calc
  IF (spdim == 1 ) THEN
      CALL calculate_moment_1d_z(vec_z,t_calc)
  ELSE IF (spdim == 2 ) THEN
      CALL calculate_moment_2d_z(vec_z,t_calc)
  ELSE IF ( spdim == 3 ) THEN
      CALL calculate_moment_3d_z(vec_z,t_calc)
  END IF
END SUBROUTINE calculate_moment_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck calculate_moment_1d_d.f
!***begin prologue     calculate_moment_1d_d
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
!***end prologue       calculate_moment_1d_d
!
  SUBROUTINE calculate_moment_1d_d(vec_1d,t_calc)
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(1))             :: vec_1d
  REAL*8                                 :: t_calc, av, avv, fac_1d
  REAL*8                                 :: exav, exavv
  INTEGER                                :: i
!
!
  av = 0.d0
  avv =0.d0
  DO i =1,nphy(1)
     fac_1d = vec_1d(i) * vec_1d(i)                           &
                        *                                     &
                 grid(1)%pt(i) 
     av = av   + fac_1d 
     avv = avv + grid(1)%pt(i) * fac_1d 
  END DO
  avv = avv - av * av
  avv = SQRT(avv)
  exav = - (beta(1) * t_calc) + x_0(1)
  exavv = .5d0 * ( sigma(1) * sigma(1) +  ( t_calc * t_calc) /  &
                                          ( sigma(1) * sigma(1) ) )
  exavv = SQRT(exavv)
  write(iout,1) t_calc, exav, av, exavv, avv
1 format(/,1x,'Time                      = ',e15.8,           &
         /,1x,'Exact <x>                 = ',e15.8,           &
         /,1x,'Calculated <x>            = ',e15.8,           &
         /,1x,'Exact <x*x> - <x><x>      = ',e15.8            &
         /,1x,'Calculated <x*x> - <x><x> = ',e15.8)
END SUBROUTINE calculate_moment_1d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck calculate_moment_2d_d.f
!***begin prologue     calculate_moment_2d_d
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
!***end prologue       calculate_moment_2d_d
!
  SUBROUTINE calculate_moment_2d_d(vec_2d,t_calc)
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(2),nphy(1))        :: vec_2d
  REAL*8, DIMENSION(2)                      :: av, avv
  REAL*8, DIMENSION(2)                      :: exav, exavv
  REAL*8                                    :: t_calc, fac_2d
  INTEGER                                   :: i, j
!
!
  av = 0.d0
  avv = 0.d0
  DO i =1,nphy(1)
     DO j=1,nphy(2)
        fac_2d = vec_2d(j,i) * vec_2d(j,i)                  &
                             *                              &
                      grid(2)%pt(j)
        av(2)  = av(2)  + fac_2d 
        avv(2) = avv(2) + fac_2d * grid(2)%pt(j)
     END DO
  END DO
  DO i =1,nphy(2)
     DO j=1,nphy(1)
        fac_2d = vec_2d(i,j) * vec_2d(i,j)                  &
                             *                              &
                     grid(1)%pt(j)
        av(1)  = av(1) + fac_2d 
        avv(1) = avv(1) + fac_2d * grid(1)%pt(j)
     END DO
  END DO
  avv(1) = avv(1) - av(1) * av(1)
  avv(2) = avv(2) - av(2) * av(2)
  avv(1) = SQRT(avv(1))
  avv(2) = SQRT(avv(2))
  exav(1) = - (beta(1) * t_calc) + x_0(1)
  exav(2) = - (beta(2) * t_calc) + x_0(2)
  exavv(1) = .5d0 * ( sigma(1) * sigma(1) + t_calc * t_calc &
            /  (sigma(1) * sigma(1)) )
  exavv(2) = .5d0 * ( sigma(2) * sigma(2) + t_calc * t_calc &
            /  (sigma(2) * sigma(2)) )
  exavv(1) = SQRT(exavv(1))
  exavv(2) = SQRT(exavv(2))
  write(iout,1) t_calc, exav(1), av(1), exavv(1), avv(1)
  write(iout,2) exav(2), av(2), exavv(2), avv(2)
1 format(/,1x,'Time                      = ',e15.8,         &
         /,1x,'Exact <x>                 = ',e15.8,         &
         /,1x,'Calculated <x>            = ',e15.8,         &
         /,1x,'Exact <x*x> - <x><x>      = ',e15.8          &
         /,1x,'Calculated <x*x> - <x><x> = ',e15.8)
2 format(/,1x,'Exact <y>                 = ',e15.8,         &
         /,1x,'Calculated <y>            = ',e15.8,         &
         /,1x,'Exact <y*y> - <y><y>      = ',e15.8          &
         /,1x,'Calculated <y*y> - <y><y> = ',e15.8)
END SUBROUTINE calculate_moment_2d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck calculate_moment_3d_d.f
!***begin prologue     calculate_moment_3d_d
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
!***end prologue       calculate_moment_3d_d
!
  SUBROUTINE calculate_moment_3d_d(vec_3d,t_calc)
  IMPLICIT NONE                          
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))     :: vec_3d
  REAL*8, DIMENSION(3)                           :: av, avv
  REAL*8, DIMENSION(3)                           :: exav, exavv
  REAL*8                                         :: fac_3d
  REAL*8                                         :: t_calc
  INTEGER                                        :: i, j, k
!
!
  av = 0.d0
  avv =0.d0
  DO i=1,nphy(1)   
     DO j =1,nphy(2)
        DO k=1,nphy(3)
           fac_3d = vec_3d(k,j,i) * vec_3d(k,j,i)                &
                                  *                              &
                           grid(3)%pt(k)
           av(3)  = av(3) + fac_3d
           avv(3) = avv(3) + fac_3d * grid(3)%pt(k)
        END DO
    END DO
  END DO
  DO i=1,nphy(1)   
     DO j =1,nphy(3)
        DO k=1,nphy(2)
           fac_3d = vec_3d(j,k,i) * vec_3d(j,k,i)                   &
                                  *                                 &
                           grid(2)%pt(k)
           av(2)  = av(2) + fac_3d
           avv(2) = avv(2) + fac_3d * grid(3)%pt(k)
        END DO
    END DO
  END DO
  DO i=1,nphy(3)   
     DO j =1,nphy(2)
        DO k=1,nphy(1)
           fac_3d = vec_3d(i,j,k) * vec_3d(i,j,k)                   &
                                  *                                 &
                           grid(1)%pt(k)
           av(1)   = av(1)  + fac_3d
           avv(1)  = avv(1) + fac_3d * grid(1)%pt(k)
        END DO
    END DO
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
  exavv(1) = .5d0 * ( sigma(1) * sigma(1) + t_calc * t_calc         &
               /  (sigma(1) * sigma(1)) )
  exavv(2) = .5d0 * ( sigma(2) * sigma(2) + t_calc * t_calc         &
               /  (sigma(2) * sigma(2)) )
  exavv(3) = .5d0 * ( sigma(3) * sigma(3) + t_calc * t_calc         &
               /  ( sigma(3) * sigma(3)) )
  exavv(1) = SQRT(exavv(1))
  exavv(2) = SQRT(exavv(2))
  exavv(3) = SQRT(exavv(3))
  write(iout,1) t_calc, exav(1), av(1), exavv(1), avv(1)
  write(iout,2) exav(2), av(2), exavv(2), avv(2)
  write(iout,3) exav(3), av(3), exavv(3), avv(3)
1 format(/,1x,'Time                      = ',e15.8,                 &
         /,1x,'Exact <x>                 = ',e15.8,                 &
         /,1x,'Calculated <x>            = ',e15.8,                 &
         /,1x,'Exact <x*x> - <x><x>      = ',e15.8                  &
         /,1x,'Calculated <x*x> - <x><x> = ',e15.8)
2 format(/,1x,'Exact <y>                 = ',e15.8,                 &
         /,1x,'Calculated <y>            = ',e15.8,                 &
         /,1x,'Exact <y*y> - <y><y>      = ',e15.8                  &
         /,1x,'Calculated <y*y> - <y><y> = ',e15.8)
3 format(/,1x,'Exact <z>                 = ',e15.8,                 &
         /,1x,'Calculated <z>            = ',e15.8,                 &
         /,1x,'Exact <z*z> - <z><z>      = ',e15.8                  &
         /,1x,'Calculated <z*z> - <z><z> = ',e15.8)
END SUBROUTINE calculate_moment_3d_d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck calculate_moment_1d_z.f
!***begin prologue     calculate_moment_1d_z
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
!***end prologue       calculate_moment_1d_z
!
  SUBROUTINE calculate_moment_1d_z(vec_1d,t_calc)
  IMPLICIT NONE                          
  COMPLEX*16, DIMENSION(nphy(1))         :: vec_1d
  COMPLEX*16                             :: conjg
  REAL*8                                 :: t_calc, av, avv, fac_1d
  REAL*8                                 :: exav, exavv
  INTEGER                                :: i
!
!
  av = 0.d0
  avv =0.d0
  DO i =1,nphy(1)
     fac_1d = vec_1d(i) * conjg( vec_1d(i) )                  &
                        *                                     &
                 grid(1)%pt(i) 
     av = av   + fac_1d 
     avv = avv + grid(1)%pt(i) * fac_1d 
  END DO
  avv = avv - av * av
  avv = SQRT(avv)
  exav = - (beta(1) * t_calc) + x_0(1)
  exavv = .5d0 * ( sigma(1) * sigma(1) +  ( t_calc * t_calc) /  &
                                          ( sigma(1) * sigma(1) ) )
  exavv = SQRT(exavv)
  write(iout,1) t_calc, exav, av, exavv, avv
1 format(/,1x,'Time                      = ',e15.8,           &
         /,1x,'Exact <x>                 = ',e15.8,           &
         /,1x,'Calculated <x>            = ',e15.8,           &
         /,1x,'Exact <x*x> - <x><x>      = ',e15.8            &
         /,1x,'Calculated <x*x> - <x><x> = ',e15.8)
END SUBROUTINE calculate_moment_1d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck calculate_moment_2d_z.f
!***begin prologue     calculate_moment_2d_z
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
!***end prologue       calculate_moment_2d_z
!
  SUBROUTINE calculate_moment_2d_z(vec_2d,t_calc)
  IMPLICIT NONE                          
  COMPLEX*16, DIMENSION(nphy(2),nphy(1))    :: vec_2d
  COMPLEX*16                                :: conjg
  REAL*8, DIMENSION(2)                      :: av, avv
  REAL*8, DIMENSION(2)                      :: exav, exavv
  REAL*8                                    :: t_calc, fac_2d
  INTEGER                                   :: i, j
!
!
  av = 0.d0
  avv = 0.d0
  DO i =1,nphy(1)
     DO j=1,nphy(2)
        fac_2d = vec_2d(j,i) * conjg ( vec_2d(j,i) )        &
                             *                              &
                      grid(2)%pt(j)
        av(2)  = av(2)  + fac_2d 
        avv(2) = avv(2) + fac_2d * grid(2)%pt(j)
     END DO
  END DO
  DO i =1,nphy(2)
     DO j=1,nphy(1)
        fac_2d = vec_2d(i,j) * conjg ( vec_2d(i,j) )        &
                             *                              &
                     grid(1)%pt(j)
        av(1)  = av(1) + fac_2d 
        avv(1) = avv(1) + fac_2d * grid(1)%pt(j)
     END DO
  END DO
  avv(1) = avv(1) - av(1) * av(1)
  avv(2) = avv(2) - av(2) * av(2)
  avv(1) = SQRT(avv(1))
  avv(2) = SQRT(avv(2))
  exav(1) = - (beta(1) * t_calc) + x_0(1)
  exav(2) = - (beta(2) * t_calc) + x_0(2)
  exavv(1) = .5d0 * ( sigma(1) * sigma(1) + t_calc * t_calc &
            /  (sigma(1) * sigma(1)) )
  exavv(2) = .5d0 * ( sigma(2) * sigma(2) + t_calc * t_calc &
            /  (sigma(2) * sigma(2)) )
  exavv(1) = SQRT(exavv(1))
  exavv(2) = SQRT(exavv(2))
  write(iout,1) t_calc, exav(1), av(1), exavv(1), avv(1)
  write(iout,2) exav(2), av(2), exavv(2), avv(2)
1 format(/,1x,'Time                      = ',e15.8,         &
         /,1x,'Exact <x>                 = ',e15.8,         &
         /,1x,'Calculated <x>            = ',e15.8,         &
         /,1x,'Exact <x*x> - <x><x>      = ',e15.8          &
         /,1x,'Calculated <x*x> - <x><x> = ',e15.8)
2 format(/,1x,'Exact <y>                 = ',e15.8,         &
         /,1x,'Calculated <y>            = ',e15.8,         &
         /,1x,'Exact <y*y> - <y><y>      = ',e15.8          &
         /,1x,'Calculated <y*y> - <y><y> = ',e15.8)
END SUBROUTINE calculate_moment_2d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck calculate_moment_3d_z.f
!***begin prologue     calculate_moment_3d_z
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
!***end prologue       calculate_moment_3d_z
!
  SUBROUTINE calculate_moment_3d_z(vec_3d,t_calc)
  IMPLICIT NONE                          
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1)) :: vec_3d
  COMPLEX*16                                     :: conjg
  REAL*8, DIMENSION(3)                           :: av, avv
  REAL*8, DIMENSION(3)                           :: exav, exavv
  REAL*8                                         :: fac_3d
  REAL*8                                         :: t_calc
  INTEGER                                        :: i, j, k
!
!
  av = 0.d0
  avv =0.d0
  DO i=1,nphy(1)   
     DO j =1,nphy(2)
        DO k=1,nphy(3)
           fac_3d = vec_3d(k,j,i) * conjg (vec_3d(k,j,i) )       &
                                  *                                &
                           grid(3)%pt(k)
           av(3)  = av(3) + fac_3d
           avv(3) = avv(3) + fac_3d * grid(3)%pt(k)
        END DO
    END DO
  END DO
  DO i=1,nphy(1)   
     DO j =1,nphy(3)
        DO k=1,nphy(2)
           fac_3d = vec_3d(j,k,i) * conjg (vec_3d(j,k,i) )          &
                                  *                                 &
                           grid(2)%pt(k)
           av(2)  = av(2) + fac_3d
           avv(2) = avv(2) + fac_3d * grid(3)%pt(k)
        END DO
    END DO
  END DO
  DO i=1,nphy(3)   
     DO j =1,nphy(2)
        DO k=1,nphy(1)
           fac_3d = vec_3d(i,j,k) * conjg ( vec_3d(i,j,k) )         &
                                  *                                 &
                           grid(1)%pt(k)
           av(1)   = av(1)  + fac_3d
           avv(1)  = avv(1) + fac_3d * grid(1)%pt(k)
        END DO
    END DO
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
  exavv(1) = .5d0 * ( sigma(1) * sigma(1) + t_calc * t_calc         &
               /  (sigma(1) * sigma(1)) )
  exavv(2) = .5d0 * ( sigma(2) * sigma(2) + t_calc * t_calc         &
               /  (sigma(2) * sigma(2)) )
  exavv(3) = .5d0 * ( sigma(3) * sigma(3) + t_calc * t_calc         &
               /  ( sigma(3) * sigma(3)) )
  exavv(1) = SQRT(exavv(1))
  exavv(2) = SQRT(exavv(2))
  exavv(3) = SQRT(exavv(3))
  write(iout,1) t_calc, exav(1), av(1), exavv(1), avv(1)
  write(iout,2) exav(2), av(2), exavv(2), avv(2)
  write(iout,3) exav(3), av(3), exavv(3), avv(3)
1 format(/,1x,'Time                      = ',e15.8,                 &
         /,1x,'Exact <x>                 = ',e15.8,                 &
         /,1x,'Calculated <x>            = ',e15.8,                 &
         /,1x,'Exact <x*x> - <x><x>      = ',e15.8                  &
         /,1x,'Calculated <x*x> - <x><x> = ',e15.8)
2 format(/,1x,'Exact <y>                 = ',e15.8,                 &
         /,1x,'Calculated <y>            = ',e15.8,                 &
         /,1x,'Exact <y*y> - <y><y>      = ',e15.8                  &
         /,1x,'Calculated <y*y> - <y><y> = ',e15.8)
3 format(/,1x,'Exact <z>                 = ',e15.8,                 &
         /,1x,'Calculated <z>            = ',e15.8,                 &
         /,1x,'Exact <z*z> - <z><z>      = ',e15.8                  &
         /,1x,'Calculated <z*z> - <z><z> = ',e15.8)
END SUBROUTINE calculate_moment_3d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE moment

