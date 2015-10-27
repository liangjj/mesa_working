!deck v_pow_exp.f90
!***begin prologue     v_pow_exp
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            r**n_p*exp(-a*r) potential
!***
!***references
!***routines called
!***end prologue       v_pow_exp
  SUBROUTINE v_pow_exp(v,pt,a,s,n_p,n,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n, n_p
  REAL*8, DIMENSION(n)                   :: v, pt
  REAL*8                                 :: a, s
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  v = v + a * pt**n_p * EXP(-s*pt)
  IF(prn) THEN
     title='a * r**n_p * exp (- s*r) potential'
     CALL prntfm(title,v,n,1,n,1,output)
  END IF
END SUBROUTINE v_pow_exp



