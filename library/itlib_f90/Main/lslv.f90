!deck lslv.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-13  Time: 16:03:58
 
!***begin prologue     lslv
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           linear system solve
!***author             schneider, barry (nsf)
!***source
!***purpose            driver for direct linear system solve.
!***
!***references

!***routines called
!***end prologue       lslv

SUBROUTINE lslv(a,b,ipvt,n,m,dim)

REAL*8, INTENT(IN OUT)                   :: a(dim,*)
REAL*8, INTENT(IN OUT)                   :: b(dim,m)
INTEGER, INTENT(IN OUT)                  :: ipvt(*)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(IN OUT)                  :: dim
IMPLICIT INTEGER (a-z)


COMMON/io/inp, iout

CALL sgetrf(n,n,a,dim,ipvt,info)
IF(info /= 0) THEN
  CALL lnkerr('error from linear solve routine')
END IF
DO  i=1,m
  CALL sgetrs('n',n,m,a,dim,ipvt,b,dim,info)
END DO
RETURN
END SUBROUTINE lslv



