!deck vexp.f
!***begin prologue     vexp
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            exponential potential
!***
!***references
!***routines called
!***end prologue       vexp
  SUBROUTINE vexp(v,pt,a,s,n,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, pt
  REAL*8                                 :: a, s
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  v = v + a*EXP(-s*pt)
  IF(prn) THEN
     title='potential'
     CALL prntfm(title,v,n,1,n,1,output)
  END IF
END SUBROUTINE vexp



