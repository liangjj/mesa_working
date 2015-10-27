! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Moment}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck moment.f
!**begin prologue     moment
!**date written       020127   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development moment
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate some moments and compare.
!**references
!**routines called
!**end prologue       moment
  SUBROUTINE moment
  USE arnoldi_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8                                 :: fpkey
  REAL*8                                 :: avx, avxx, exavx, exavxx
  REAL*8, DIMENSION(3)                   :: new_sigma
  CHARACTER (LEN=2)                      :: itoc
  LOGICAL                                :: dollar
  CHARACTER (LEN=80)                     :: chrkey
  INTEGER                                :: i
  IF(spdim /= 1) THEN
     WRITE(iout,1)
     CALL lnkerr('quit')
  END IF
  IF( dollar('$initial-state',card,title,inp) ) THEN
     CALL fparr(card,'alpha',alpha,3,' ')
     CALL fparr(card,'sigma',sigma,3,' ')
     CALL fparr(card,'x_0',x_0,3,' ')
     CALL fparr(card,'beta',beta,3,' ')
  END IF
  new_sigma = sigma/sqrt(alpha)
  chi = chi + psi0
  avx=0.d0
  avxx=0.d0
  IF(typke == 'fd') THEN
     DO  i=1,nphy(1)
         avx = avx + grid(1)%wt(i) * grid(1)%pt(i) * chi(i) * CONJG(chi(i))
         avxx = avxx + grid(1)%wt(i) * grid(1)%pt(i) * &
                       grid(1)%pt(i) * chi(i) * CONJG(chi(i))
     END DO
  ELSE
     DO  i=1,nphy(1)
         avx = avx   + grid(1)%pt(i) * chi(i) * CONJG(chi(i))
         avxx = avxx + grid(1)%pt(i) * grid(1)%pt(i) * chi(i) * CONJG(chi(i))
     END DO
  END IF
  avxx = avxx - avx*avx
  avxx=SQRT(avxx)
  exavx=-(beta(1)*t1) + x_0(1)
  exavxx = .5d0 * ( new_sigma(1) * new_sigma(1) + t1 * t1 &
               /  (new_sigma(1)*new_sigma(1)) )
  exavxx=SQRT(exavxx)
  WRITE(iout,2) exavx, avx, exavxx, avxx
  chi = chi - psi0
1    FORMAT(/,1X,'This is a test for one-dimensional free wavepacket',  &
    /,1X,'Current calculation does not qualify')
2    FORMAT(/,1X,'Exact <x>                 = ',e15.8,  &
    /,1X,'Calculated <x>            = ',e15.8,  &
    /,1X,'Exact <x*x> - <x><x>      = ',e15.8,  &
    /,1X,'Calculated <x*x> - <x><x> = ',e15.8)
END SUBROUTINE moment
