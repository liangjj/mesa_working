!deck lavec.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2008-03-01  Time: 16:31:57
 
!***begin prologue     lavec
!***date written       970531   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           davidson, trial, vector
!***author             schneider, barry (nsf)
!***source
!***purpose            set up new trial vector based on some zeroth
!***                   order model.  this is equivalent to a preconditioning
!***                   of the matrix.
!***references

!***routines called
!***end prologue       lavec

SUBROUTINE lavec(energy,diag,resid,n,m,iter,prnt)

REAL*8, INTENT(IN)                       :: energy
REAL*8, INTENT(IN)                       :: diag(n)
REAL*8, INTENT(OUT)                      :: resid(n,m)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(IN OUT)                  :: iter
LOGICAL, INTENT(IN)                      :: prnt
IMPLICIT INTEGER (a-z)

REAL*8 test, zero, nrzero, one

CHARACTER (LEN=4) :: itoc
CHARACTER (LEN=80) :: title
DATA zero, nrzero, one / 0.d0, 1.0D-06, 1.d0 /

COMMON/io/inp, iout

DO  i=1,m
  DO  j=1,n
    test = energy - diag(j)
    IF(ABS(test) >= nrzero) THEN
      resid(j,i) = resid(j,i)/test
    ELSE
      resid(j,i) = one
    END IF
  END DO
END DO
IF(prnt) THEN
  title='new trial vectors iteration = '//itoc(iter)
  CALL prntrm(title,resid,n,m,n,m,iout)
END IF
RETURN
END SUBROUTINE lavec





