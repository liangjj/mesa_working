!deck maket.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:45:11
 
!***begin prologue     maket
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            symmetrized kinetic energy matrix elements.
!***
!***description        computes the matrix of second derivatives
!                      plus block operator.
!***references

!***routines called
!***end prologue       maket

SUBROUTINE maket(fa,dfb,ddfb,pt,wt,secder,bloch,mat,typarg,n)

REAL*8, INTENT(IN)                       :: fa(n,n)
REAL*8, INTENT(IN OUT)                   :: dfb(n,n)
REAL*8, INTENT(IN OUT)                   :: ddfb(n,n)
REAL*8, INTENT(IN OUT)                   :: pt(n)
REAL*8, INTENT(IN)                       :: wt(n)
REAL*8, INTENT(IN OUT)                   :: secder(n,n)
REAL*8, INTENT(IN OUT)                   :: bloch(n,n)
REAL*8, INTENT(OUT)                      :: mat(n,n)
CHARACTER (LEN=*), INTENT(IN)            :: typarg
INTEGER, INTENT(IN)                      :: n
IMPLICIT INTEGER (a-z)




COMMON/io/inp, iout

IF(typarg == 'linear') THEN
  CALL copy(secder,mat,n*n)
  CALL vadd(mat,mat,bloch,n*n)
ELSE IF(typarg == 'quadratic') THEN
  DO  i=1,n
    DO  j=1,n
      mat(i,j) = fa(i,i) * wt(i) * 4.d0 * ( pt(i)*pt(i)*ddfb(i,j) + dfb(i,j) )
    END DO
  END DO
  CALL vadd(mat,mat,bloch,n*n)
ELSE
  CALL lnkerr('error in argument type')
END IF
RETURN
END SUBROUTINE maket



