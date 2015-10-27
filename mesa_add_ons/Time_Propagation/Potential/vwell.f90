!deck vwell.f
!***begin prologue     vwell
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            potential well
!***
!***references
!***routines called
!***end prologue       vwell
  SUBROUTINE vwell(v,depth,n,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v
  REAL*8                                 :: depth
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  INTEGER                                :: i
  DO  i=1,n
      v(i) = v(i) + depth
  END DO
  IF(prn) THEN
     title='potential'
     CALL prntrm(title,v,n,1,n,1,iout)
  END IF
END SUBROUTINE vwell



