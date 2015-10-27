!deck addang.f
!***begin prologue     addang
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            centrifugal potential
!***
!***references
!***routines called
!***end prologue       addang

  SUBROUTINE addang(v,pt,scale,l,coord,n,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: v, pt
  REAL*8                                 :: scale
  INTEGER                                :: l, lnum
  CHARACTER (LEN=*)                      :: coord
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  lnum=l*(l+1)
  IF(coord == 'rho') THEN
     lnum=l*l
  END IF
  v = v - scale*lnum/(pt*pt)
  IF(prn) THEN
     title='potential l = '//itoc(l)
     CALL prntfm(title,v,n,1,n,1,iout)
   END IF
END SUBROUTINE addang



