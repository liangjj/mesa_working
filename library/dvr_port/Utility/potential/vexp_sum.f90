!deck vexp_sum.f
!***begin prologue     vexp_sum
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            sum of two exponentials
!***
!***references
!***routines called
!***end prologue       vexp_sum
  SUBROUTINE vexp_sum(v,pt,a,s,n,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, pt
  REAL*8, DIMENSION(2)                   :: a, s
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  v = v + a(1) * EXP(-s(1)*pt*pt) + a(2)* EXP(-s(2)*pt*pt) 
  IF(prn) THEN
     title='sum of two exponentials potential'
     CALL prntfm(title,v,n,1,n,1,output)
  END IF
END SUBROUTINE vexp_sum



