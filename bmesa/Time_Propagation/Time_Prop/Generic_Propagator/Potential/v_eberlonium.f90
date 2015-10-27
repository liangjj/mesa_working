!deck v_eberlonium.f
!***begin prologue     v_eberlonium
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            eberlonium potential
!***
!***references
!***routines called
!***end prologue       v_eberlonium
  SUBROUTINE v_eberlonium(v,pt,z,a,b,npwr,n,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n, npwr
  REAL*8, DIMENSION(n)                   :: v, pt
  REAL*8                                 :: z, a, b
  LOGICAL                                :: prn
  CHARACTER (len=80)                     :: title
  v = v + z / sqrt( (a*(pt**npwr) + b) )
  IF(prn) THEN
     title='potential'
     CALL prntrm(title,v,n,1,n,1,iout)
  END IF
END SUBROUTINE v_eberlonium



