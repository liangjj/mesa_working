! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Trial}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck trial
!**begin prologue     trial
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Arnoldi, trial
!***
!**author             schneider, b. i.(nsf)
!**source             trial
!**purpose            trial vectors for  Arnoldi method.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       trial
  SUBROUTINE trial
  USE arnoldi_global
  IMPLICIT NONE
  INTEGER                                :: i, nout 
!
!  First vector is taken as solution from previous step.  The space is
!  completed with unit vectors and the entire set schmidt orthonormalized.
!
  vec(:,1) = psi0
  IF(ntrial > 1) THEN
     vec(:,2:ntrial)=(0.d0,0.d0)
      DO  i=2,ntrial
          vec(i,i)=(1.d0,0.d0)
      END DO
  END IF
  CALL cschmt(vec,thresh,n3d,1,ntrial,nout,.true.,.false.)
END SUBROUTINE trial






















