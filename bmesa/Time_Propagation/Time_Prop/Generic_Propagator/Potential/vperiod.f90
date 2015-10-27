!deck vperiod.f
!***begin prologue     vperiod
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            periodic potential
!***
!***references
!***routines called
!***end prologue       vperiod
  SUBROUTINE vperiod(v,pt,a,n,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, pt
  REAL*8                                 :: a
  LOGICAL                                :: prn
  CHARACTER (len=80)                     :: title
  v = v + a * ( 1.d0 - cos(pt) )
  IF(prn) THEN
     title='periodic potential'
     CALL prntrm(title,v,n,1,n,1,iout)
  END IF
END SUBROUTINE vperiod



