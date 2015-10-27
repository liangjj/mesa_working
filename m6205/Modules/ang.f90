!deck ang.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2014-06-29  Time: 11:13:38
 
!***begin prologue     ang
!***date written       930802   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           ang, link m6200
!***author             schneider, barry (nsf)
!***source             m6200
!***purpose            fill angular quantities in lebedev quadrature
!***references

!***routines called
!***end prologue       ang

SUBROUTINE ang (pt,cthet,sthet,cphi,sphi,phpt,nleb)

REAL*8, INTENT(IN)                       :: pt(3,nleb)
REAL*8, INTENT(OUT)                      :: cthet(nleb)
REAL*8, INTENT(OUT)                      :: sthet(nleb)
REAL*8, INTENT(OUT)                      :: cphi(nleb)
REAL*8, INTENT(OUT)                      :: sphi(nleb)
REAL*8, INTENT(OUT)                      :: phpt(nleb)
INTEGER, INTENT(IN)                      :: nleb
IMPLICIT INTEGER (a-z)
REAL*8  twopi, rsq, r


COMMON /io/ inp, iout
DATA twopi /  6.283185307179586D0  /

DO  i=1,nleb
  rsq=pt(1,i)*pt(1,i)+pt(2,i)*pt(2,i)+pt(3,i)*pt(3,i)
  r=SQRT(rsq)
  cthet(i)=pt(3,i)/r
  sthet(i)=SQRT(1.d0-cthet(i)*cthet(i))
  IF (sthet(i) /= 0.d0) THEN
    cphi(i)=pt(1,i)/(r*sthet(i))
    sphi(i)=pt(2,i)/(r*sthet(i))
  ELSE
    cphi(i)=1.d0
    sphi(i)=0.d0
  END IF
  phpt(i)=ATAN2(sphi(i),cphi(i))
END DO
RETURN
END SUBROUTINE ang

