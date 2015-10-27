!deck cheb2.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:40:31
 
!***begin prologue     cheb2
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           chebyshev weights
!***author             schneider, b. i.(nsf)
!***source
!***purpose            weight functions and their first and second
!***                   derivatives
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       cheb2

SUBROUTINE cheb2(wt,dwt,ddwt,x,n)


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
  fac = ( 1.d0 - x(i) * x(i) ) ** .25D0
  wt(i) = fac
  inv2 = inv * inv
  dwt(i) = - .5D0 * x(i) * fac * inv
  ddwt(i) = -.5D0 * fac * inv * x(i) *x(i) * fac *inv2
END DO
RETURN
END SUBROUTINE cheb2















