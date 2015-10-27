!deck cpy_3.f
!***begin prologue     cpy_3
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           tridiagonal matrix fill
!***author             schneider, b. i.(nsf)
!***source
!***purpose            fill a tridiagonal matrix
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       cpy_3
  SUBROUTINE cpy_3(band,e,d,n)
  IMPLICIT NONE
  INTEGER                   :: n
  REAL*8, DIMENSION(2,n)    :: band
  REAL*8, DIMENSION(n)      :: e, d
  e = band(1,:)
  d = band(2,:)
END SUBROUTINE cpy_3
