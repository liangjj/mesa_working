!deck lsolve.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2008-03-01  Time: 16:32:29
 
!***begin prologue     lsolve
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           linear system solve
!***author             schneider, barry (nsf)
!***source
!***purpose            driver for direct linear system solve.
!***
!***references

!***routines called
!***end prologue       lsolve

SUBROUTINE lsolve(a,b,ipvt,n,m,dim)

REAL*8, INTENT(IN OUT)                   :: a(dim,*)
REAL*8, INTENT(IN OUT)                   :: b(dim,m)
INTEGER, INTENT(IN OUT)                  :: ipvt(*)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(IN OUT)                  :: dim
IMPLICIT INTEGER (a-z)


COMMON/io/inp, iout

CALL sgefa(a,dim,n,ipvt,info)
IF(info /= 0) THEN
  CALL lnkerr('error from linear solve routine')
END IF
DO  i=1,m
  CALL sgesl(a,dim,n,ipvt,b(1,i),0)
END DO
RETURN
END SUBROUTINE lsolve
