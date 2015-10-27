!deck vwel12.f
!***begin prologue     vwel12
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            two dimensional potential well
!***
!***references
!***routines called
!***end prologue       vwel12
  SUBROUTINE vwel12(v,q1,q2,a12,d12,n1,n2,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n1, n2 
  REAL*8, DIMENSION(n2,n1)               :: v
  REAL*8, DIMENSION(n1)                  :: q1
  REAL*8, DIMENSION(n2)                  :: q2
  REAL*8                                 :: a12, d12
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  INTEGER                                :: i, j
  DO  i=1,n1
      DO  j=1,n2
          IF( ABS(q1(i)-q2(j)) <= a12) THEN
              v(j,i) = v(j,i) + d12
          END IF
      END DO
  END DO
  IF(prn) THEN
     title='two dimensional well'
     CALL prntrm(title,v,n2,n1,n2,n1,iout)
  END IF
1    FORMAT(/,9X,'coordinate',8X,'potential')
2    FORMAT(5X,e15.8,3X,e15.8)
END SUBROUTINE vwel12



