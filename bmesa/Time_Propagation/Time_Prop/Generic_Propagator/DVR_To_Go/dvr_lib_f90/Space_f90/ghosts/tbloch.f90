!deck tbloch.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:46:04
 
!***begin prologue     tbloch
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            lobatto matrix elements of bloch
!***                   operator.
!***
!***description        computes the matrix

!                      M = - fa(xright) * fb'(xright)
!                          + fa(xleft)  * fb'(xleft)
!***references

!***routines called
!***end prologue       tbloch

SUBROUTINE tbloch(fa,dfb,pt,mat,typarg,n)

REAL*8, INTENT(IN)                       :: fa(n,n)
REAL*8, INTENT(IN)                       :: dfb(n,n)
REAL*8, INTENT(IN)                       :: pt(n)
REAL*8, INTENT(OUT)                      :: mat(n,n)
CHARACTER (LEN=*), INTENT(IN)            :: typarg
INTEGER, INTENT(IN)                      :: n
IMPLICIT INTEGER (a-z)



COMMON/io/inp, iout

CALL rzero(mat,n*n)
IF(typarg == 'linear') THEN
  DO  i=1,n
    mat(n,i) = mat(n,i) - fa(n,n)*dfb(n,i)
    mat(1,i) = mat(1,i) + fa(1,1)*dfb(1,i)
  END DO
ELSE IF(typarg == 'quadratic') THEN
  DO  i=1,n
    mat(n,i) = mat(n,i) - 2.d0*pt(n)*pt(n)*fa(n,n)*dfb(n,i)
    mat(i,i) = mat(1,i) + 2.d0*pt(1)*pt(1)*fa(1,1)*dfb(1,i)
  END DO
ELSE
  CALL lnkerr('error in argument type')
END IF
RETURN
END SUBROUTINE tbloch



