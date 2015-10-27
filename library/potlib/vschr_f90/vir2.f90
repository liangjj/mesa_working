!deck vir2.f
!***begin prologue     vir2
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            inverse square
!***
!***references
!***routines called
!***end prologue       vir2
  SUBROUTINE vir2(v,pt,n,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, pt
  LOGICAL                                :: prn
  REAL*8                                 :: one=1.d0
  CHARACTER (LEN=80)                     :: title
  v = v + one/( pt * pt )
  IF(prn) THEN
     title='potential'
     CALL prntfm(title,v,n,1,n,1,iout)
  END IF
END SUBROUTINE vir2


