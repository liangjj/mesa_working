!deck vir4.f
!***begin prologue     vir4
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            potential well
!***
!***references
!***routines called
!***end prologue       vir4
  SUBROUTINE vir4(v,pt,n,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, pt
  LOGICAL                                :: prn
  REAL*8                                 :: one=1.d0
  CHARACTER (LEN=80)                     :: title
  v = v + one/( one + pt )**4
  IF(prn) THEN
     title='potential'
     CALL prntfm(title,v,n,1,n,1,output)
  END IF
END SUBROUTINE vir4



