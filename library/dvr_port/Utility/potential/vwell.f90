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
  SUBROUTINE vwell(v,pt,LEN,depth,n,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, pt
  REAL*8                                 :: len, depth
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  INTEGER                                :: i
  DO  i=1,n
      IF(( pt(i) - pt(1) ) <= LEN) THEN
         v(i) = v(i) + depth
      END IF
  END DO
  IF(prn) THEN
     title='potential'
     CALL prntrm(title,v,n,1,n,1,output)
  END IF
END SUBROUTINE vwell



