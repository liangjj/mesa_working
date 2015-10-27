!deck newsol.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-13  Time: 16:04:08
 
!***begin prologue     newsol
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           best current solution
!***author             schneider, barry (nsf)
!***source
!***purpose            express best current davidson solution in original
!***                   basis.
!***
!***references

!***routines called
!***end prologue       newsol

SUBROUTINE newsol(sol,vec,coef,n,m,nrhs,maxvec)

REAL*8, INTENT(IN OUT)                   :: sol(n,nrhs)
REAL*8, INTENT(IN OUT)                   :: vec(n,maxvec)
REAL*8, INTENT(IN OUT)                   :: coef(maxvec,nrhs)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: m
INTEGER, INTENT(IN OUT)                  :: nrhs
INTEGER, INTENT(IN OUT)                  :: maxvec
IMPLICIT INTEGER (a-z)


COMMON/io/inp, iout

CALL ebcxx(sol,vec,coef,n,m,nrhs,n,n,maxvec)
RETURN
END SUBROUTINE newsol



