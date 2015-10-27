!deck jacobi.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:43:55
 
!***begin prologue     jacobi
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           jacobi weights
!***author             schneider, b. i.(nsf)
!***source
!***purpose            weight functions and their first and second
!***                   derivatives
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       jacobi

SUBROUTINE jacobi(wt,dwt,ddwt,x,alf,bet,n)


REAL*8, INTENT(OUT)                      :: wt(n)
REAL*8, INTENT(OUT)                      :: dwt(n)
REAL*8, INTENT(OUT)                      :: ddwt(n)
REAL*8, INTENT(IN)                       :: x(n)
REAL*8, INTENT(IN)                       :: alf
REAL*8, INTENT(IN)                       :: bet
INTEGER, INTENT(IN)                      :: n
IMPLICIT INTEGER (a-z)
REAL*8  inv, inv1, fac, fac1

COMMON/io/inp, iout

DO  i=1,n
  inv = 1.d0 / ( 1.d0 - x(i) )
  inv1 = 1.d0 / ( 1.d0 + x(i) )
  fac = ( 1.d0 - x(i) )**(.5*alf )
  fac1 = ( 1.d0 + x(i) )**(.5*bet )
  wt(i) = fac * fac1
  dwt(i) = - .5D0 * alf * fac * inv * fac1 + fac * .5D0 * bet * fac1 * inv1
  ddwt(i) = -.5D0 * alf * ( .5D0 * alf - 1.d0 ) * fac * fac1 *  &
      inv * inv + .5D0 * bet * ( .5D0 * bet - 1.d0 ) * fac * fac1 *inv1 * inv1
END DO
RETURN
END SUBROUTINE jacobi















