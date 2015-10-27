!deck h12h0.f
!***begin prologue     h12h0
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            transform a 1-D vector from the dvr to
!***                   the H0 representation.
!***
!***references

!***routines called
!***end prologue       h12h0

  SUBROUTINE h12h0(vecout,vecin)
  USE dvr_shared
  USE dvd_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),addvec)     :: vecin, vecout
  CALL ebtc(vecout,grid(1)%eigvec,vecin,nphy(1),nphy(1),addvec)
END SUBROUTINE h12h0
