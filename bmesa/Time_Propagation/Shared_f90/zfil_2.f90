! \documentclass{article}
! Code converted using TO_F90 by Alan Miller
! Date: 2003-02-06  Time: 07:16:42
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Zfil_2}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck zfil_2.f
!**begin prologue     zfil_2
!**date written       960723   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            fill array
!**references
!**routines called
!**end prologue       zfil_2
  SUBROUTINE zfil_2(z1,z2)
  USE arnoldi_global
  USE dvr_shared
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(nphy(1))         :: z1
  COMPLEX*16, DIMENSION(nphy(2))         :: z2
  INTEGER                                :: i, j, count 
  count=0
  DO  i=1,nphy(1)
      DO  j=1,nphy(2)
          count=count+1
          psi0(count) = z1(i)*z2(j)
      END DO
  END DO
END SUBROUTINE zfil_2


