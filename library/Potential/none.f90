!deck none.f
!***begin prologue     none
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            no potential
!***
!***references
!***routines called
!***end prologue       none
  SUBROUTINE NONE(v,pt,n,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, pt
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  IF(prn) THEN
     title='potential'
     CALL prntrm(title,v,n,1,n,1,iout)
  END IF
END SUBROUTINE NONE



