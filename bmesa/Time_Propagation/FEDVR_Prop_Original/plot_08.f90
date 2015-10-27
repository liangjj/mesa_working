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
!deck plot_08
!***begin prologue     plot_08
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            combine initial solution with solution to
!                      driven equation to get total solution.
!***
!***references
!***routines called
!***end prologue       plot_08
!
  SUBROUTINE plot_08(t)
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: t, i
  CHARACTER (LEN=16)                     :: fptoc
  REAL*8                                 :: rtemp
!
!
  IF (log_main(8)) THEN
      title='solution at t = '//fptoc(tim_pts(t+1))
      CALL plot_wavefunction(psi_08)
  END IF
  IF (plot.and.t == ntreg) THEN
      IF(spdim==1) then
         CALL plot_1d(psi_08)
      ELSE IF(spdim==2) THEN
         CALL plot_2d(psi_08)
      ELSE
         CALL plot_3d(psi_08)
      END IF
  END IF
  call iosys('read real "initial state" from bec',2*n3d,v_scr_08,0,' ')
  IF(spdim==1) THEN
     CALL auto_08_1d(auto_corr,v_scr_08,psi_08)
  ELSE IF(spdim==2) THEN
     CALL auto_08_2d(auto_corr,v_scr_08,psi_08,fac)
  ELSE IF(spdim==3) THEN
     CALL auto_08_3d(auto_corr,v_scr_08,psi_08,fac,v_1)
  END IF
END SUBROUTINE plot_08
