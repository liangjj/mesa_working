!deck h1e.f
!***begin prologue     h1e
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            energy scale a 1-D vector.
!***
!***references

!***routines called
!***end prologue       h1e

  SUBROUTINE h1e(vecout)
  USE dvr_shared
  USE dvd_global
  IMPLICIT NONE
  INTEGER                                :: i, j
  REAL*8                                 :: test 
  REAL*8, DIMENSION(nphy(1),addvec)      :: vecout
  DO  i=1,addvec
      DO  j=1,nphy(1)
          test = eig(i) - grid(1)%eigv(j)
          IF(ABS(test) >= nrzero) THEN
             vecout(j,i) = vecout(j,i)/test
          ELSE
             vecout(j,i) = one
          END IF
      END DO
  END DO
END SUBROUTINE h1e
