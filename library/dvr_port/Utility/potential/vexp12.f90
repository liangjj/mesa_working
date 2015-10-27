!deck vexp12.f
!***begin prologue     vexp12
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            exponential pair potential
!***
!***references
!***routines called
!***end prologue       vexp12
  SUBROUTINE vexp12(v,q1,q2,a12,b12,n1,n2,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n1, n2
  REAL*8, DIMENSION(n2,n1)               :: v
  REAL*8, DIMENSION(n1)                  :: q1
  REAL*8, DIMENSION(n2)                  :: q2
  REAL*8                                 :: a12, b12
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  INTEGER                                :: i
  DO  i=1,n1
      v(:,i) = v(:,i) + a12 * EXP ( - b12 * ABS( q2(:) - q1(i) ) )
  END DO
  IF(prn) THEN
     title='two dimensional exponential interaction'
     CALL prntrm(title,v,n2,n1,n2,n1,output)
  END IF
END SUBROUTINE vexp12



