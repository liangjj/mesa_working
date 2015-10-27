!deck dfadfb.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:41:27
 
!***begin prologue     dfadfb
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            lobatto matrix elements of first derivative
!***                   of function with first derivative of function.
!***
!***references

!***routines called
!***end prologue       dfadfb

SUBROUTINE dfadfb(dfa,dfb,pt,wt,mat,typarg,n)

REAL*8, INTENT(IN)                       :: dfa(n,n)
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
      DO  k=1,n
        mat(i,j) = mat(i,j) + dfa(k,i)*wt(k)*dfb(k,j)
      END DO
    END DO
  END DO
ELSE IF(typarg == 'quadratic') THEN
  DO  i=1,n
    DO  j=1,n
      DO  k=1,n
        mat(i,j) = mat(i,j) + 4.d0*pt(k)*pt(k)*wt(k)* dfa(k,i)*dfb(k,j)
      END DO
    END DO
  END DO
ELSE
  CALL lnkerr('error in argument type')
END IF
RETURN
END SUBROUTINE dfadfb









