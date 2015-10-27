!deck vcoul.f
!***begin prologue     vcoul
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            potential well
!***
!***references
!***routines called
!***end prologue       vcoul
  SUBROUTINE vcoul(v,pt,z,n,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, pt
  REAL*8                                 :: z
  LOGICAL                                :: prn
  CHARACTER (len=80)                     :: title
  v = v + z / pt
  IF(prn) THEN
     title='potential'
     CALL prntrm(title,v,n,1,n,1,output)
  END IF
END SUBROUTINE vcoul



