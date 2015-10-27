!deck h12dvr.f
!***begin prologue     h12dvr
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
!***end prologue       h12dvr
  SUBROUTINE h12dvr(vecout,vecin)
  USE dvr_shared
  USE dvd_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),addvec)      :: vecout, vecin
  CALL ebc(vecout,grid(1)%eigvec,vecin,nphy(1),nphy(1),addvec)
END SUBROUTINE h12dvr
