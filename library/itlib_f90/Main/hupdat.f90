!deck hupdat.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-13  Time: 16:03:30
 
!***begin prologue     hupdat
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose
!***
!***references

!***routines called
!***end prologue       hupdat

SUBROUTINE hupdat(h,htmp,b,btmp,vec,hvec,rhs,n,nrhs,nbeg,nend,m)

REAL*8, INTENT(OUT)                      :: h(m,m)
REAL*8, INTENT(OUT)                      :: htmp(m,m)
REAL*8, INTENT(OUT)                      :: b(m,*)
REAL*8, INTENT(OUT)                      :: btmp(m,*)
REAL*8, INTENT(IN)                       :: vec(n,*)
REAL*8, INTENT(IN)                       :: hvec(n,*)
REAL*8, INTENT(IN)                       :: rhs(n,nrhs)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: nrhs
INTEGER, INTENT(IN)                      :: nbeg
INTEGER, INTENT(IN)                      :: nend
INTEGER, INTENT(IN OUT)                  :: m
IMPLICIT INTEGER (a-z)
REAL*8  sdot


COMMON/io/inp, iout

DO  i=1,nend
  DO  j=nbeg,nend
    h(i,j) = sdot(n,vec(1,i),1,hvec(1,j),1)
    h(j,i) = sdot(n,vec(1,j),1,hvec(1,i),1)
  END DO
  DO  j=1,nrhs
    b(i,j) = sdot(n,vec(1,i),1,rhs(1,j),1)
  END DO
END DO
DO  i=1,nend
  DO  j=1,i
    htmp(i,j) = h(i,j)
    htmp(j,i) = h(j,i)
  END DO
  DO  j=1,nrhs
    btmp(i,j) = b(i,j)
  END DO
END DO
RETURN
END SUBROUTINE hupdat

