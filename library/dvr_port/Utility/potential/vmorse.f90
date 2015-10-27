!deck vmorse.f
!***begin prologue     vmorse
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            morse potential
!***
!***references
!***routines called
!***end prologue       vmorse
  SUBROUTINE vmorse(v,pt,a,b,exp_1,exp_2,e_0,n,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, pt
  REAL*8                                 :: a, b, exp_1, exp_2, e_0
  LOGICAL                                :: prn
  CHARACTER (len=80)                     :: title
  v = v + 4.d0 * a * ( exp(-exp_1*pt) - 2.d0 * exp(exp_2*pt) ) + b &
        - e_0 * pt
  IF(prn) THEN
     title='potential'
     CALL prntrm(title,v,n,1,n,1,output)
  END IF
END SUBROUTINE vmorse



