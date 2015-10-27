! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{C_vect0}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck c_vect0.f
!**begin prologue     c_vect0
!**date written       960615   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           hamiltonian
!**author             schneider, barry (nsf)
!**source             3-dim
!**purpose            projection of zeroth order eigenvectors on to
!**                   dvr basis.
!**references
!**routines called
!**end prologue       c_vect0
  SUBROUTINE c_vect0(etemp,itemp,root)
  USE arnoldi_global
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                                 :: root
  REAL*8, DIMENSION(n3d)                  :: etemp
  INTEGER, DIMENSION(n3d,*)               :: itemp
  REAL*8                                  :: tmp
  CHARACTER (LEN=16)                      :: fptoc
  INTEGER                                 :: i, j, k
  INTEGER                                 :: count
!
!       The eig and ind arrays are destroyed by this routine.
!
  psi0=(0.d0,0.d0)
  IF(spdim == 1) THEN
     energy=grid(1)%eigv_0(root)
     psi0 = grid(1)%eigvec_0(:,root)
  ELSE IF(spdim == 2) THEN
     call psi0_fil_2(etemp,itemp,root)     
     energy=etemp(root)
  ELSE IF(spdim == 3) then
     call psi0_fil_3(etemp,itemp,root)     
     energy=etemp(root)
  END IF
  WRITE(iout,1) root, energy
1    FORMAT(/,5X,'initial state = ',i3,/,5X, 'energy        = ',e15.8)
END SUBROUTINE c_vect0






