!deck lagure.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:44:29
 
!***begin prologue     lagure
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           laguree weights
!***author             schneider, b. i.(nsf)
!***source
!***purpose            weight functions and their first and second
!***                   derivatives
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       lagure

SUBROUTINE lagure(wt,dwt,ddwt,x,alf,n)


REAL*8, INTENT(OUT)                      :: wt(n)
REAL*8, INTENT(OUT)                      :: dwt(n)
REAL*8, INTENT(OUT)                      :: ddwt(n)
REAL*8, INTENT(IN)                       :: x(n)
REAL*8, INTENT(IN)                       :: alf
INTEGER, INTENT(IN)                      :: n
IMPLICIT INTEGER (a-z)
REAL*8  fac, inv, inv2, exfac

COMMON/io/inp, iout

DO  i=1,n
  fac = x(i)**(.5D0*alf)
  inv = 1.d0 / x(i)
  inv2 = inv * inv
  exfac = EXP(-.5D0*x(i))
  wt(i) = fac*exfac
  dwt(i) = .5D0 * alf * fac * inv * exfac - .5D0 * fac * exfac
  ddwt(i) = ( .5D0 * alf * ( .5D0 * alf - 1.d0 ) * fac * inv2  &
      - .5D0 * alf * fac * inv + .25D0 * fac ) * exfac
END DO
RETURN
END SUBROUTINE lagure















