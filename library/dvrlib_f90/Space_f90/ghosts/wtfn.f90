!deck wtfn.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:46:18
 
!***begin prologue     wtfn
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           wtfn functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            weight functions and their first and second
!***                   derivatives
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       wtfn

SUBROUTINE wtfn(px,pwt,TYPE,alpha,beta,ngot,n)


INTEGER*8, INTENT(IN OUT)                :: px
INTEGER*8, INTENT(IN OUT)                :: pwt
CHARACTER (LEN=*), INTENT(IN)            :: TYPE
REAL*8, INTENT(IN OUT)                   :: alpha
REAL*8, INTENT(IN OUT)                   :: beta
INTEGER, INTENT(IN OUT)                  :: ngot
INTEGER, INTENT(IN)                      :: n
IMPLICIT INTEGER (a-z)

REAL*8 x, wt
#ifdef decpointer

#END IF decpointer
#ifdef sgipointer

#END IF sgipointer
COMMON/io/inp, iout
pointer (px,x(1))
pointer (pwt,wt(1))

need=wptoin(3*n)
CALL memory(need,pwt,ngot,'weights',0)
p2wt=1
p2dwt=p2wt+n
p2ddwt=p2dwt+n
IF(TYPE == 'legendre') THEN
  CALL vfill(wt(p2wt),1.d0,n)
  CALL rzero(wt(p2dwt),n)
  CALL rzero(wt(p2ddwt),n)
ELSE IF(TYPE == 'hermite') THEN
  CALL hermit(wt(p2wt),wt(p2dwt),wt(p2ddwt),x,n)
ELSE IF(TYPE == 'laguerre') THEN
  CALL lagure(wt(p2wt),wt(p2dwt),wt(p2ddwt),x,alpha,n)
ELSE IF(TYPE == 'chebyshev-1') THEN
  CALL cheb1(wt(p2wt),wt(p2dwt),wt(p2ddwt),x,n)
ELSE IF(TYPE == 'chebyshev-2') THEN
  CALL cheb2(wt(p2wt),wt(p2dwt),wt(p2ddwt),x,n)
ELSE IF(TYPE == 'jacobi') THEN
  CALL cheb2(wt(p2wt),wt(p2dwt),wt(p2ddwt),x,alpha,beta,n)
ELSE
  CALL lnkerr('quit. error in weight type')
END IF
RETURN
END SUBROUTINE wtfn
















