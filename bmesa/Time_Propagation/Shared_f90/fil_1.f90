!deck fil_1.f
!***begin prologue     fil_1
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose
!***references
!***routines called
!***end prologue       fil_1
  SUBROUTINE fil_1(phi,n)
  USE arnoldi_global 
  IMPLICIT NONE
  INTEGER                                :: n
  REAL*8, DIMENSION(n)                   :: phi
  psi0=phi
END SUBROUTINE fil_1

