!deck fil_08.f
!***begin prologue     fil_08
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose
!***references
!***routines called
!***end prologue       fil_08
  SUBROUTINE fil_08(phix,phiy,phiz,nx,ny,nz)
  USE dvrprop_global
  IMPLICIT NONE
  INTEGER                                :: nx, ny,nz
  REAL*8, DIMENSION(nx)                  :: phix
  REAL*8, DIMENSION(ny)                  :: phiy
  REAL*8, DIMENSION(nz)                  :: phiz
  INTEGER                                :: i, j, k, count
  count = 0
  DO  i=1,nx
    DO  j=1,ny
        do k=1,nz
           count = count + 1
           psi0(count) = phix(i) * phiy(j) * phiz(k)
        END DO
    END DO
  END DO
END SUBROUTINE fil_08