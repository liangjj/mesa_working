!deck vanhmo.f
!***begin prologue     vanhmo
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            anharmonic oscillator
!***
!***references
!***routines called
!***end prologue       vhmo
  SUBROUTINE vanhmo(v,pt,n,prn)
  USE inout
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, pt
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  v = v + .5D0*( pt*pt + pt*pt*pt*pt )
  IF(prn) THEN
     title='potential'
     CALL prntfm(title,v,n,1,n,1,output)
  END IF
END SUBROUTINE vanhmo



