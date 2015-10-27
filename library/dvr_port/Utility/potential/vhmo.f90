!deck vhmo.f
!***begin prologue     vhmo
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            potential well
!***
!***references
!***routines called
!***end prologue       vhmo
  SUBROUTINE vhmo(v,pt,fac,n,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n 
  REAL*8, DIMENSION(n)                   :: v, pt
  REAL*8                                 :: fac
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  v = v + fac*pt*pt
  IF(prn) THEN
     title='potential'
     CALL prntfm(title,v,n,1,n,1,output)
  END IF
END SUBROUTINE vhmo



