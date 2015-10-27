!deck h2e.f
!***begin prologue     h2e
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            energy scale a 2-D vector.
!***
!***references

!***routines called
!***end prologue       h2e
  SUBROUTINE h2e(vecout)
  USE dvd_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                                     :: i, j, k
  REAl*8, DIMENSION(nphy(2),nphy(1),addvec)   :: vecout
  REAL*8  fac
  REAL*8 test
!
  DO  i=1,addvec
      DO  j=1,nphy(1)
          DO  k=1,nphy(2)
              test = eig(i) - grid(1)%eigv(j) - grid(2)%eigv(k)
              IF(ABS(test) >= nrzero) THEN
                 vecout(k,j,i) = vecout(k,j,i)/test
              ELSE
                 vecout(k,j,i) = one
              END IF
          END DO
      END DO
  END DO
END SUBROUTINE h2e
