!deck rwell.f
!***begin prologue     rwell
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            rounded potential well
!***
!***references
!***routines called
!***end prologue       rwell
  SUBROUTINE rwell(v,pt,nwell,awell,n,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, pt
  INTEGER                                :: nwell
  REAL*8                                 :: awell
  LOGICAL                                :: prn
  REAL*8                                 :: one=1.d0
  CHARACTER (LEN=80)                     :: title
  v = v - one/SQRT( awell + pt**nwell )
  IF(prn) THEN
     title='potential for this region'
     CALL prntrm(title,v,n,1,n,1,iout)
  END IF
END SUBROUTINE rwell



