!deck fafb.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:42:45
 
!***begin prologue     fafb
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            lobatto overlap matrix elements.
!***
!***references

!***routines called
!***end prologue       fafb

SUBROUTINE fafb(fa,fb,wt,mat,n)

REAL*8, INTENT(IN)                       :: fa(n,n)
REAL*8, INTENT(IN)                       :: fb(n,n)
REAL*8, INTENT(IN)                       :: wt(n)
REAL*8, INTENT(OUT)                      :: mat(n,n)
INTEGER, INTENT(IN)                      :: n
IMPLICIT INTEGER (a-z)


COMMON/io/inp, iout

CALL rzero(mat,n*n)
DO  i=1,n
  mat(i,i) = fa(i,i)*wt(i)*fb(i,i)
END DO
RETURN
END SUBROUTINE fafb
