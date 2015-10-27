!deck h3e.f
!***begin prologue     h3e
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            energy scale a 2-D vector.
!***
!***references

!***routines called
!***end prologue       h3e
  SUBROUTINE h3e(vecout)
  USE dvd_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                                             :: i, j, k, l
  REAl*8, DIMENSION(nphy(3),nphy(2),nphy(1),addvec)   :: vecout
  REAL*8                                             :: test 
  DO  i=1,addvec
      DO  j=1,nphy(1)
          DO  k=1,nphy(2)
              DO  l=1,nphy(3)
                  test = eig(i) - grid(1)%eigv(j)    &
                                - grid(2)%eigv(k)    &
                                - grid(3)%eigv(l)
                  IF(ABS(test) >= nrzero) THEN
                     vecout(l,k,j,i) = vecout(l,k,j,i)/test
                  ELSE
                     vecout(l,k,j,i) = one
                  END IF
              END DO
          END DO
      END DO
  END DO
RETURN
END SUBROUTINE h3e
