! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Subroutine for Radial Wavepacket}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck rad_paket.f
!**begin prologue     rad_paket
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate zero time rad wavepacket
!**references
!**routines called
!**end prologue       rad_paket
  SUBROUTINE rad_paket
  USE io
  USE prop_global
  USE arnoldi_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  COMPLEX*16                             :: norm, cdotc
  INTEGER                                :: i  
  WRITE(iout,1)
  IF(typ_pak == 'exponential') then
     WRITE(iout,2)
     psi0 = ( grid(1)%pt - x_0(1) )**powr(1) * &
              exp( - ( sigma(1) + eye*beta(1) ) * ( grid(1)%pt - x_0(1) ) )
  ELSE IF(typ_pak == 'gaussian') then
     WRITE(iout,3)
     psi0 = ( grid(1)%pt - x_0(1) )**powr(1) * &
            exp( - sigma(1) * ( grid(1)%pt - x_0(1) ) * &
                              ( grid(1)%pt - x_0(1) )   &
                 - eye * ( grid(1)%pt - x_0(1) ) )
  END IF
  WRITE(iout,4)
  WRITE(iout,5) sigma(1), x_0(1), beta(1), powr(1)
  psi0 = psi0 * grid(1)%wt
  do i=1,nphy(1)
     psi0(i) = psi0(i) * grid(1)%f(i,i)
  END DO
  norm = 1.d0/sqrt(cdotc(nphy(1),psi0,1,psi0,1))
  psi0 = psi0 * norm  
  norm = cdotc(nphy(1),psi0,1,psi0,1)
  title='initial wavepacket'
  write(iout,7) norm
  call prntcm(title,psi0,nphy(1),1,nphy(1),1,iout)
1    FORMAT(/,5X,'initial wavepacket at t=0')
2    FORMAT(/,5X,'the form of the radial packet is:',///,5X,  &
                 'psi = (r - x_0 )**n * exp( - ( alpha + i*beta ) * ' &
                                           '( r - x_0 ))')
3    FORMAT(/,5X,'the form of the radial packet is:',///,5X,  &
                 'psi = ( r - x_0 )**n * &
                  exp( - alpha * ( r -x_0 ) * ( r - x_0 ) - i * beta * ' &
                                                            '(r- x_0 ) )')
4    FORMAT(/,1X,'radial wave packet parameters')
5    FORMAT(/,1X,'exponent    = ',e15.8, &
            /,1X,'shift       = ',e15.8, &
            /,1X,'momentum    = ',e15.8, &
            /,1X,'power       = ',i3)
7    FORMAT(/,2x,'Normalization of Initial Wavepacket = ',e15.8,1x,e15.8)
END SUBROUTINE rad_paket







