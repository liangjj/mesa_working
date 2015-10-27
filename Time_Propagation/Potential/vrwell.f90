!deck vrwell.f
!***begin prologue     vrwell
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            potential well
!***
!***references
!***routines called
!***end prologue       vrwell
  SUBROUTINE vrwell(v,pt,nwell,awell,n,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n, nwell
  REAL*8, DIMENSION(n)                   :: v, pt
  REAL*8                                 :: awell
  LOGICAL                                :: prn
  REAL*8                                 :: spt, one=1.d0
  CHARACTER (LEN=80)                     :: title
  v = v - one/SQRT( awell + pt**nwell )
  IF(prn) THEN
     title='potential'
     CALL prntfm(title,v,n,1,n,1,iout)
  END IF
END SUBROUTINE vrwell



