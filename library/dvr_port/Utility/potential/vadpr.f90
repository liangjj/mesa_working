!deck vadpr.f
!***begin prologue     vadpr
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            adiabatic pair potential
!***
!***references
!***routines called
!***end prologue       vadpr
  SUBROUTINE vadpr(v,q1,q2,a12,b12,n,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, q1, q2
  REAL*8                                 :: a12, b12
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  v = v + a12 * EXP ( - b12 * ABS( q2 - q1 ) )
  IF(prn) THEN
     title='two dimensional exponential interaction'
     CALL prntrm(title,v,n,1,n,1,output)
  END IF
END SUBROUTINE vadpr



