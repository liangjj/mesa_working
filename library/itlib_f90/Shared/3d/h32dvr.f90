!deck h32dvr.f
!***begin prologue     h32dvr
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            transform a 3-D vector from the dvr to
!***                   the H0 representation.
!***
!***references

!***routines called
!***end prologue       h32dvr
  SUBROUTINE h32dvr(vecout,vecin,tmp,tmp1)
  USE dvd_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                                             :: i, j    
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),addvec)   :: vecin
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),addvec)   :: vecout
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),addvec)   :: tmp
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),addvec)   :: tmp1
  CALL ebc(tmp1,grid(3)%eigvec,vecin,                          &
           nphy(3),nphy(3),nphy(2)*nphy(1)*addvec)
  DO  i=1,nphy(1)
      DO  j=1,addvec
          CALL ebct(tmp(1,1,i,j),tmp1(1,1,i,j),grid(2)%eigvec, &
                    nphy(3),nphy(2),nphy(2))
      END DO
  END DO
  DO  i=1,addvec
      CALL ebct(vecout(1,1,1,i),tmp(1,1,1,i),grid(1)%eigvec,   &
                nphy(3)*nphy(2),nphy(1),nphy(1))
  END DO
END SUBROUTINE h32dvr
