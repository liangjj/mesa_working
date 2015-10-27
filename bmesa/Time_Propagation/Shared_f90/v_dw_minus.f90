!deck v_dw_minus.f90
!***begin prologue     v_dw_minus
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            
!***
!***references
!***routines called
!***end prologue       v_dw_minus
  SUBROUTINE v_dw_minus(v,q,t,n,prn)
  USE io
  USE potential
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v
  REAL*8, DIMENSION(n)                   :: q
  REAL*8                                 :: t, x_dw
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  INTEGER                                :: i
  x_dw = (x_final - x_initial) * ( t - t_initial)/ (t_final - t_initial) + x_initial
  v = v - nu_dw_minus * EXP ( - ( q(:) - x_dw) * ( q - x_dw ) / ( gamma_dw * gamma_dw ) )
  IF(prn) THEN
     title='two dimensional exponential interaction'
     CALL prntrm(title,v,n,n,n,n,iout)
  END IF
END SUBROUTINE v_dw_minus



