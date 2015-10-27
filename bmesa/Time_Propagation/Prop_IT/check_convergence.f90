!***********************************************************************
                           MODULE check_convergence
                           INTERFACE chk_con
                    MODULE PROCEDURE chk_con_d, chk_con_z
                       END INTERFACE chk_con
!***********************************************************************
                           CONTAINS
!***********************************************************************
!***********************************************************************
!***********************************************************************
!**deck chk_con_z
!**begin prologue     chk_con_z
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           propagation, arnoldi
!**author             schneider, barry (nsf)
!**source
!**purpose            Module to check convergence
!**description        Check the difference between the initial and
!**                   computed vector.  Finished is issued if converged.
!**references
!**routines called
!***end prologue       chk_con_z
  SUBROUTINE chk_con_z(v_out)
  USE arnoldi_global_rt
  IMPLICIT NONE
  COMPLEX*16,DIMENSION(n3d)              :: v_out
  COMPLEX*16                             :: cdotc
  REAL*8                                 :: norm
  INTEGER                                :: i
  norm = 0.d0
  vscr = v_out - vscr
  norm = sqrt ( cdotc(n3d,vscr,1,vscr,1) )
  WRITE(iout,1) norm
  cntrl='continue'
  IF(norm <= cnverg) THEN
     cntrl='finished'
  ELSE
     vscr = v_out
  END IF
1    FORMAT(/,1X,'overlap               = ',e15.8)
END SUBROUTINE chk_con_z
!**deck chk_con_d
!**begin prologue     chk_con_d
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           propagation, arnoldi
!**author             schneider, barry (nsf)
!**source
!**purpose            convergence check
!**description
!**references
!**routines called
!***end prologue       chk_con_d
  SUBROUTINE chk_con_d(v_out)
  USE arnoldi_global_it
  IMPLICIT NONE
  REAL*8,DIMENSION(n3d)              :: v_out
  REAL*8                             :: ddot
  REAL*8                             :: norm
  INTEGER                            :: i
  norm = 0.d0
  vscr = v_out - vscr
  norm = sqrt ( ddot(n3d,vscr,1,vscr,1) )
  WRITE(iout,1) norm
  cntrl='continue'
  IF(norm <= cnverg) THEN
     cntrl='finished'
  ELSE
     vscr = v_out
  END IF
1    FORMAT(/,1X,'overlap               = ',e15.8)
END SUBROUTINE chk_con_d
END MODULE check_convergence
