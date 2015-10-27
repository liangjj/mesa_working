!deck fill_d.f
!***begin prologue     fill_d
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose
!***references
!***routines called
!***end prologue       fill_d
  SUBROUTINE fill_d(phix,phiy,phiz,nx,ny,nz)
  USE dvrprop_global_rt
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
           wave_function(count,1) = phix(i) * phiy(j) * phiz(k)
        END DO
    END DO
  END DO
END SUBROUTINE fill_d
