! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Diagonalize 3 Point FD Matrix Elements}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck moment_08.f
!***begin prologue     moment_08
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
!***end prologue       moment_08
!
  SUBROUTINE moment_08
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE                          
  REAL*8                                 :: fpkey
  REAL*8, DIMENSION(3)                   :: av, avv, exav, exavv
  REAL*8, DIMENSION(3)                   :: new_sigma
  CHARACTER (LEN=2)                      :: itoc
  LOGICAL                                :: dollar
  CHARACTER (LEN=80)                     :: chrkey
  INTEGER                                :: i
!
!
  IF ( dollar('$initial-state',card,title,inp) ) THEN
       CALL fparr(card,'alpha',alpha,3,' ')
       CALL fparr(card,'sigma',sigma,3,' ')
       CALL fparr(card,'x_0',x_0,3,' ')
       CALL fparr(card,'beta',beta,3,' ')
  END IF
  new_sigma = sigma/sqrt(alpha)
  IF(spdim == 1) THEN
     call mom_1d(psi_08,av,avv)
     avv(1) = avv(1) - av(1) * av(1)
     avv(1) = SQRT(avv(1))
     exav(1) = - (beta(1)*t1) + x_0(1)
     exavv(1) = .5d0 * ( new_sigma(1) * new_sigma(1) + t1 * t1 &
               /  (new_sigma(1)*new_sigma(1)) )
     exavv(1) = SQRT(exavv(1))
     write(iout,2) exav(1), av(1), exavv(1), avv(1)
  ELSE IF(spdim == 2) THEN
     call mom_2d(psi_08,fac,av,avv)
     avv(1) = avv(1) - av(1) * av(1)
     avv(2) = avv(2) - av(2) * av(2)
     avv(1) = avv(1) - av(1) * av(1)
     avv(2) = avv(2) - av(2) * av(2)
     avv(1) = SQRT(avv(1))
     avv(2) = SQRT(avv(2))
     exav(1) = - (beta(1)*t1) + x_0(1)
     exav(2) = - (beta(2)*t1) + x_0(2)
     exavv(1) = .5d0 * ( new_sigma(1) * new_sigma(1) + t1 * t1 &
               /  (new_sigma(1)*new_sigma(1)) )
     exavv(2) = .5d0 * ( new_sigma(2) * new_sigma(2) + t1 * t1 &
               /  (new_sigma(2)*new_sigma(2)) )
     exavv(1) = SQRT(exavv(1))
     exavv(2) = SQRT(exavv(2))
     write(iout,2) exav(1), av(1), exavv(1), avv(1)
     write(iout,3) exav(2), av(2), exavv(2), avv(2)
  ELSE IF(spdim == 3) THEN
     call mom_3d(psi_08,fac,v_1,av,avv)
     avv(1) = avv(1) - av(1) * av(1)
     avv(2) = avv(2) - av(2) * av(2)
     avv(3) = avv(3) - av(3) * av(3)
     avv = abs(avv)
     avv(1) = SQRT(avv(1))
     avv(2) = SQRT(avv(2))
     avv(3) = SQRT(avv(3))
     exav(1) = - (beta(1)*t1) + x_0(1)
     exav(2) = - (beta(2)*t1) + x_0(2)
     exav(3) = - (beta(3)*t1) + x_0(3)
     exavv(1) = .5d0 * ( new_sigma(1) * new_sigma(1) + t1 * t1 &
               /  (new_sigma(1)*new_sigma(1)) )
     exavv(2) = .5d0 * ( new_sigma(2) * new_sigma(2) + t1 * t1 &
               /  (new_sigma(2)*new_sigma(2)) )
     exavv(3) = .5d0 * ( new_sigma(3) * new_sigma(3) + t1 * t1 &
               /  (new_sigma(3)*new_sigma(3)) )
     exavv(1) = SQRT(exavv(1))
     exavv(2) = SQRT(exavv(2))
     exavv(3) = SQRT(exavv(3))
     write(iout,2) exav(1), av(1), exavv(1), avv(1)
     write(iout,3) exav(2), av(2), exavv(2), avv(2)
     write(iout,4) exav(3), av(3), exavv(3), avv(3)
  END IF
1 format(/,1x,'This is a test for a separable free wavepacket', &
         /,1x,'Current calculation does not qualify')
2 format(/,1x,'Exact <x>                 = ',e15.8,      &
         /,1x,'Calculated <x>            = ',e15.8,      &
         /,1x,'Exact <x*x> - <x><x>      = ',e15.8       &
         /,1x,'Calculated <x*x> - <x><x> = ',e15.8)
3 format(/,1x,'Exact <y>                 = ',e15.8,      &
         /,1x,'Calculated <y>            = ',e15.8,      &
         /,1x,'Exact <y*y> - <y><y>      = ',e15.8       &
         /,1x,'Calculated <y*y> - <y><y> = ',e15.8)
4 format(/,1x,'Exact <y>                 = ',e15.8,      &
         /,1x,'Calculated <y>            = ',e15.8,      &
         /,1x,'Exact <z*z> - <z><z>      = ',e15.8       &
         /,1x,'Calculated <z*z> - <z><z> = ',e15.8)
END SUBROUTINE moment_08
