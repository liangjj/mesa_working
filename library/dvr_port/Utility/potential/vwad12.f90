!deck vwad12.f
!***begin prologue     vwad12
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            two dimensional potential well in adiabatic
!***                   coordinates
!***
!***references
!***routines called
!***end prologue       vwad12
  SUBROUTINE vwad12(v,q1,q2,depth,n,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, q1, q2
  REAL*8                                 :: depth
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  WRITE(output,1) depth
  v = v + depth
  IF(prn) THEN
     title='two dimensional well'
     CALL prntrm(title,v,n,1,n,1,output)
  END IF
1    FORMAT(/,1X,'two dimensional potential well depth = ',e15.8)
END SUBROUTINE vwad12



