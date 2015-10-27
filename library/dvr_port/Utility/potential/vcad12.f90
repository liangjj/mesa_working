!deck vcad12.f
!***begin prologue     vcad12
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            one-d model coulomb interaction
!***
!***references
!***routines called
!***end prologue       vcad12
  SUBROUTINE vcad12(v,x,y,z,n,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, x, y
  REAL*8                                 :: z
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  WRITE(output,1) z
  v = v + z / (x + y)
  IF(prn) THEN
     title='potential'
     CALL prntrm(title,v,n,1,n,1,output)
  END IF
1    FORMAT(/,1X,'potential = z / ( x1 + x2 ) where',/,1X,  &
    '         z = ',e15.8)
END SUBROUTINE vcad12



