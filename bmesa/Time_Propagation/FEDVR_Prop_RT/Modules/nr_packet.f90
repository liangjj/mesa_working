! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Subroutine for Gaussian Wavepacket Normalization}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck nr_packet.f
!**begin prologue     nr_packet
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            calculate normailzation integrals for gaussian wavepacket
!**references
!**routines called
!**end prologue       nr_packet
  SUBROUTINE nr_packet(norm)
  USE arnoldi_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8                                 :: norm
  REAL*8, DIMENSION(3)                   :: ov 
  INTEGER                                :: i
  IF(system == 'cartesian') THEN
     norm=1.d0
     DO i=1,spdim
        call ov1_quad(ov(i),grid(i)%pt,grid(i)%wt, &
                      x_0(i),alpha(i),sigma(i),nphy(i))
        norm=norm*ov(i)
     END DO
  ELSE
     call ov2_quad(norm,grid(1)%pt,grid(1)%wt, &
                   x_0(1),alpha(1),sigma(1),nphy(1))    
  END IF
  norm=1.d0/SQRT(norm)
  WRITE(iout,1) norm
1    FORMAT(/,1X,'overlap integral for gaussian packet = ',e15.8)
END SUBROUTINE nr_packet






