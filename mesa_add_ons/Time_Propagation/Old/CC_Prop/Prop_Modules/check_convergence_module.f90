!***********************************************************************
                           MODULE check_convergence_module
                           USE Iterative_Global
!***********************************************************************
                           INTERFACE check_conv
                    MODULE PROCEDURE chk_con_d, chk_con_z
                       END INTERFACE check_conv
!***********************************************************************
                           CONTAINS
!***********************************************************************
!***********************************************************************
!***********************************************************************
!**deck chk_con_d
!**begin prologue     chk_con_d
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           propagation, lanczos
!**author             schneider, barry (nsf)
!**source
!**purpose            convergence check
!**description
!**references
!**routines called
!***end prologue       chk_con_d
  SUBROUTINE chk_con_d(v_out)
  IMPLICIT NONE
  REAL*8,DIMENSION(:)                :: v_out
  REAL*8                             :: ddot
  REAL*8                             :: norm
  INTEGER                            :: i
  norm = 0.d0
  vscr_d = v_out - vscr_d
  norm = sqrt ( ddot(n3d,vscr_d,1,vscr_d,1) )
  WRITE(iout,1) norm
  cntrl='continue'
  IF(norm <= cnverg) THEN
     cntrl='finished'
  ELSE
     vscr_d = v_out
  END IF
1    FORMAT(/,1X,'overlap               = ',e15.8)
END SUBROUTINE chk_con_d
!***********************************************************************
!***********************************************************************
!**deck chk_con_z
!**begin prologue     chk_con_z
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           propagation, lanczos
!**author             schneider, barry (nsf)
!**source
!**purpose            Module to check convergence
!**description        Check the difference between the initial and
!**                   computed vector.  Finished is issued if converged.
!**references
!**routines called
!***end prologue       chk_con_z
  SUBROUTINE chk_con_z(v_out)
  IMPLICIT NONE
  COMPLEX*16,DIMENSION(:)                :: v_out
  COMPLEX*16                             :: cdotc
  REAL*8                                 :: norm
  INTEGER                                :: i
  norm = 0.d0
  vscr_z = v_out - vscr_z
  norm = sqrt ( cdotc(n3d,vscr_z,1,vscr_z,1) )
  WRITE(iout,1) norm
  cntrl='continue'
  IF(norm <= cnverg) THEN
     cntrl='finished'
  ELSE
     vscr_z = v_out
  END IF
1    FORMAT(/,1X,'overlap               = ',e15.8)
END SUBROUTINE chk_con_z
!***********************************************************************
!***********************************************************************
END MODULE check_convergence_module
