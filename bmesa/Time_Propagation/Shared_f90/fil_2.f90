!deck fil_2.f
!***begin prologue     fil_2
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose
!***references
!***routines called
!***end prologue       fil_2
  SUBROUTINE fil_2(phix,phiy,nx,ny)
  USE arnoldi_global
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: nx, ny
  REAL*8, DIMENSION(nx)                  :: phix
  REAL*8, DIMENSION(ny)                  :: phiy
  INTEGER                                :: i, j, count
  count = 0
  DO  i=1,nx
      DO j=1,ny
         count = count + 1
         psi0(count) = phix(i) * phiy(j)
      END DO
  END DO
END SUBROUTINE fil_2

