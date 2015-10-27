!deck h22dvr.f
!***begin prologue     h22dvr
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            transform a 2-D vector from the dvr to
!***                   the H0 representation.
!***
!***references

!***routines called
!***end prologue       h22dvr

  SUBROUTINE h22dvr(vecout,vecin,tmp)
  USE dvd_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                                     :: i    
  REAL*8, DIMENSION(nphy(2),nphy(1),addvec)   :: vecin
  REAL*8, DIMENSION(nphy(2),nphy(1),addvec)   :: vecout
  REAL*8, DIMENSION(nphy(2),nphy(1),addvec)   :: tmp
  CALL ebc(tmp,grid(2)%eigvec,vecin,                       &
           nphy(2),nphy(2),nphy(1)*addvec)
  DO  i=1,addvec
      CALL ebct(vecout(1,1,i),tmp(1,1,i),grid(1)%eigvec,   &
                nphy(2),nphy(1),nphy(1))
  END DO
END SUBROUTINE h22dvr
