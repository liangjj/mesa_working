!deck vonv.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-13  Time: 16:04:44
 
!***begin prologue     vonv
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***routines called
!***end prologue       vonv

SUBROUTINE vonv(vecout,v,vecin,n,nc,m)

REAL*8, INTENT(OUT)                      :: vecout(n,nc,m)
REAL*8, INTENT(IN)                       :: v(n,nc,nc)
REAL*8, INTENT(IN)                       :: vecin(n,nc,m)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: nc
INTEGER, INTENT(IN)                      :: m
IMPLICIT INTEGER (a-z)
REAL*8  one

COMMON/io/inp, iout

CALL rzero(vecout,n*nc*m)
DO  i=1,n
  DO  j=1,m
    DO  k=1,nc
      DO  l=1,nc
        vecout(i,k,j) = vecout(i,k,j) + v(i,k,l)*vecin(i,l,j)
      END DO
    END DO
  END DO
END DO
RETURN
END SUBROUTINE vonv
