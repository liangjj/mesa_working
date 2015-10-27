!deck cheb1.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:40:15
 
!***begin prologue     cheb1
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           cheb1 weights
!***author             schneider, b. i.(nsf)
!***source
!***purpose            weight functions and their first and second
!***                   derivatives
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       cheb1

SUBROUTINE cheb1(wt,dwt,ddwt,x,n)


REAL*8, INTENT(OUT)                      :: wt(n)
REAL*8, INTENT(OUT)                      :: dwt(n)
REAL*8, INTENT(OUT)                      :: ddwt(n)
REAL*8, INTENT(IN)                       :: x(n)
INTEGER, INTENT(IN)                      :: n
IMPLICIT INTEGER (a-z)
REAL*8  inv, inv2

COMMON/io/inp, iout

DO  i=1,n
  inv = 1.d0 / ( 1.d0 - x(i) * x(i) )
  fac = inv ** .25D0
  inv2 = inv * inv
  wt(i) = fac
  dwt(i) = .5D0 * x(i) * fac * inv
  ddwt(i) = .5D0 * fac * inv + 1.25D0 * x(i) * x(i) * fac * inv2
END DO
RETURN
END SUBROUTINE cheb1















