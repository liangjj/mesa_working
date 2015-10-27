!deck fadfb.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:42:33
 
!***begin prologue     fadfb
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            lobatto matrix elements of function
!***                   and first derivative.
!***
!***references

!***routines called
!***end prologue       fadfb

SUBROUTINE fadfb(fa,dfb,pt,wt,mat,typarg,n)

REAL*8, INTENT(IN)                       :: fa(n,n)
REAL*8, INTENT(IN)                       :: dfb(n,n)
REAL*8, INTENT(IN)                       :: pt(n)
REAL*8, INTENT(IN)                       :: wt(n)
REAL*8, INTENT(OUT)                      :: mat(n,n)
CHARACTER (LEN=*), INTENT(IN)            :: typarg
INTEGER, INTENT(IN)                      :: n
IMPLICIT INTEGER (a-z)



COMMON/io/inp, iout

CALL rzero(mat,n*n)
IF(typarg == 'linear') THEN
  DO  i=1,n
    DO  j=1,n
      mat(i,j) = fa(i,i)*wt(i)*dfb(i,j)
    END DO
  END DO
ELSE IF(typarg == 'quadratic') THEN
  DO  i=1,n
    DO  j=1,n
      mat(i,j) = fa(i,i)*wt(i)*dfb(i,j)*2.d0*pt(i)
    END DO
  END DO
ELSE
  CALL lnkerr('error in argument type')
END IF
RETURN
END SUBROUTINE fadfb



