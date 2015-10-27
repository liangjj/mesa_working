!***********************************************************************
                           MODULE exponential_propagator
                           INTERFACE exp_prop
                    MODULE PROCEDURE exp_prop_d, exp_prop_z
                       END INTERFACE exp_prop
!***********************************************************************
                           CONTAINS
!***********************************************************************
!***********************************************************************
!***********************************************************************
!deck exp_prop_z.f
!**begin prologue     exp_prop_z
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           propagation
!**author             schneider, barry (nsf)
!**source
!**purpose            exponentiation of hamiltonian matrix
!**
!**description
!**references
!**routines called
!**end prologue       exp_prop_z
 SUBROUTINE exp_prop_z(v_in,v_out)
  USE arnoldi_global_rt,           v_int => work
  IMPLICIT NONE
  INTEGER                                :: i, j
  COMPLEX*16, DIMENSION(n3d)             :: v_in, v_out
  COMPLEX*16                             :: cdotc
  REAL*8                                 :: norm
  LOGICAL, DIMENSION(4)                  :: prnt
  WRITE(iout,1) deltat
!
  CALL cehbtc(v_int,u,v_in,size,n3d,1)
  v_int(1:size) = EXP(-eye*eig(1:size)*deltat/hbar) * v_int(1:size)
  CALL cebc(v_out,u,v_int,n3d,size,1)
  IF(drct == 'subtract') THEN
     v_out = v_out - v_in
  END IF
1    FORMAT(/,1X,'exponentiating at delta t = ',e15.8)
END SUBROUTINE exp_prop_z
!deck exp_prop_d.f
!**begin prologue     exp_prop_d
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           propagation
!**author             schneider, barry (nsf)
!**source
!**purpose            exponentiation of hamiltonian matrix
!**
!**description
!**references
!**routines called
!**end prologue       exp_prop_d
 SUBROUTINE exp_prop_d(v_in,v_out)
  USE arnoldi_global_it,           v_int => work
  IMPLICIT NONE
  INTEGER                                :: i, j
  REAL*8, DIMENSION(n3d)                 :: v_in, v_out
  REAL*8                                 :: ddot
  REAL*8                                 :: norm
  LOGICAL, DIMENSION(4)                  :: prnt
  WRITE(iout,1) deltat
!
  CALL ebtc(v_int,u,v_in,size,n3d,1)
  v_int(1:size) = EXP(-eig(1:size)*deltat/hbar) * v_int(1:size)
  CALL ebc(v_out,u,v_int,n3d,size,1)
  norm= 1.d0/SQRT( ddot(n3d,v_out,1,v_out,1) )
  v_out = norm * v_out
  IF(drct == 'subtract') THEN
     v_out = v_out - v_in
  END IF
1    FORMAT(/,1X,'exponentiating at delta t = ',e15.8)
END SUBROUTINE exp_prop_d
END MODULE exponential_propagator
