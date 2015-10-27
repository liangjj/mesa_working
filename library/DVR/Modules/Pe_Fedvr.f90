!deck Pe_Fedvr.f
!***begin prologue     Pe_Fedvr
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate the parts of the regional matrices
!***                   not depending on the m quantum number.  since a different
!***                   form for DVR basis functions are needed for even and odd
!***                   m, two different kinetic energy operators need to be constructed.
!***                   Later, the even and odd m kinetic energy will be built from these.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Pe_Fedvr

  SUBROUTINE Pe_Fedvr(grid,type)
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  TYPE (coordinates)        :: grid
  CHARACTER(LEN=*)          :: type
  REAL*8                    :: dum
  INTEGER                   :: i
!
!
  IF ( type == 'eta') THEN
       DO i = nreg
          ALLOCATE( grid%reg_vec(i)%vr( npt(i) ) )
          grid%reg_vec(i)%vr(:) = R_ab * ( z_b - z_a )                     &
                                       *                                   &
                                    grid%reg_pt_wt(i)%qr(:)
       END DO
!
  ELSE IF ( type == 'xi' ) THEN
       DO i = nreg
          ALLOCATE( grid%reg_vec(i)%vr( npt(i) ) )
          grid%reg_vec(i)%vr(:) = R_ab * ( z_a + z_b )                     &
                                       *                                   &
                                  grid%reg_pt_wt(i)%qr(:)
       END DO
  END IF
!
END SUBROUTINE Pe_Fedvr
