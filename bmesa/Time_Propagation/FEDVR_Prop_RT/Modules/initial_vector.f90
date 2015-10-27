!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE initial_vector
!**begin prologue     initial_vector
!**date written       960615   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           hamiltonian
!**author             schneider, barry (nsf)
!**source             3-dim
!**purpose            Sort eigenvectors according energy and place
!***                  the desired root into initial state.
!**references
!**routines called
!**end prologue       initial_vector
  INTERFACE c_vect
  SUBROUTINE c_vect_2d(wave_function,etemp,itemp,root)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                 :: root
  REAL*8, DIMENSION(nphy(2),nphy(1))      :: wave_function
  REAL*8, DIMENSION(nphy(2)*nphy(1))      :: etemp
  INTEGER, DIMENSION(nphy(2)*nphy(1),2)   :: itemp
  REAL*8                                  :: tmp
  INTEGER                                 :: i, j, i1, j1
  INTEGER                                 :: ii, count
  END SUBROUTINE c_vect_2d
  SUBROUTINE c_vect_3d(wave_function,etemp,itemp,root)
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                         :: root
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1))      :: wave_function
  REAL*8, DIMENSION(nphy(3)*nphy(2)*nphy(1))      :: etemp
  INTEGER, DIMENSION(nphy(3)*nphy(2)*nphy(1),3)   :: itemp
  REAL*8                                          :: tmp
  INTEGER                                         :: i, j, k, i1, j1, k1
  INTEGER                                         :: ii, count
  END SUBROUTINE c_vect_3d
  END INTERFACE c_vect
END MODULE initial_vector

