!deck v_spheroidal.f
!***begin prologue     v_spheroidal
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            exponential potential
!***
!***references
!***routines called
!***end prologue       v_spheroidal
  SUBROUTINE v_spheroidal(v,pt,f,Z_a,Z_b,R_ab,n,coord,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: pt
  REAL*8, DIMENSION(n)                   :: v
  REAL*8, DIMENSION(n,n)                 :: f
  REAL*8                                 :: Z_a
  REAL*8                                 :: Z_b
  REAL*8                                 :: R_ab
  LOGICAL                                :: prn
  INTEGER                                :: i
  CHARACTER (LEN=*)                      :: coord
  CHARACTER (LEN=80)                     :: title
  IF (coord == 'eta' ) THEN
      DO i = 1, n
         v(i) = v(i) + pt(i) * R_ab * ( Z_b - Z_a ) * f(i,i) * f(i,i)
      END DO
  ELSE IF (coord == 'xi' ) THEN
      DO i = 1, n
         v(i) = v(i) + pt(i) * R_ab * ( Z_a + Z_b ) * f(i,i) * f(i,i)
      END DO
  END IF
  IF(prn) THEN
     title='potential'
     CALL prntfm(title,v,n,1,n,1,iout)
  END IF
END SUBROUTINE v_spheroidal



