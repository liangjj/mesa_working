!deck vres.f
!***begin prologue     vres
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            resonant exponential potential
!***
!***references
!***end prologue       vres
  SUBROUTINE vres(v,pt,a,s,shift,n,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                      :: n
  REAL*8, DIMENSION(n)         :: v, pt
  REAL*8, DIMENSION(2)         :: a, s
  REAL*8                       :: shift
  LOGICAL                      :: prn
  CHARACTER (LEN=80)           :: title
  v = v + a(1)*EXP(-s(1)*pt) +  &
          a(2)*EXP(-s(2)*( pt - shift )*( pt - shift ) )
  IF(prn) THEN
     title='potential'
     CALL prntfm(title,v,n,1,n,1,iout)
  END IF
END SUBROUTINE vres



