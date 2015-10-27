!deck Re_KE.f
!***begin prologue     Re_KE
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            1. calculate piecewise lobatto dvr functions and
!***                      their one-body matrices
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Re_KE

  SUBROUTINE Re_KE(tr, left_end_mat_el, right_end_mat_el, inv_sqrt_wt, n , region)
  USE dvr_global,            ONLY : nreg
  IMPLICIT NONE
  REAL*8, DIMENSION(:,:)          :: tr
  REAL*8                          :: end_mat_el
  REAL*8, DIMENSION (:)           :: inv_sqrt_wt
  INTEGER                         :: n
  INTEGER                         :: region
  INTEGER                         :: i
  INTEGER                         :: j
!
!
!
  IF ( nreg == 1 ) THEN
!
!      No bridge functions, so its straightforward.
!
       DO i = 1, n
          DO j = 1, i
             tr(i,j) = inv_sqrt_wt(i) * tr(i,j) * inv_sqrt_wt(j)
             tr(j,i) = tr(i,j)
          END DO
       END DO
       return
  END IF
!
  IF ( region == 1 ) THEN
!
!      There is a bridge function at the right hand boundary that requires the
!      the matrix element involving the first function in the next sector.  This
!      is in the right_end_mat_el.
!      
       DO i = 1, n - 1
          DO j = 1, i
             tr(i,j) = inv_sqrt_wt(i) * tr(i,j) * inv_sqrt_wt(j)
             tr(j,i) = tr(i,j)
          END DO
       END DO
       i = n
       DO j = 1, n - 1
          tr(i,j) = inv_sqrt_wt(i) * tr(i,j) * inv_sqrt_wt(j)
          tr(j,i) = tr(i,j)
       END DO
       tr(i,i) = inv_sqrt_wt(i) * ( tr(i,i) + right_end_mat_el ) * inv_sqrt_wt(i)
!
  ELSE IF ( region == nreg )
!
!      This is the last sector.  Here we require information involving the bridge function
!      at the left hand boundary.
!
       i = 1
       tr(i,i) = inv_sqrt_wt(i) * ( tr(i,i) + left_end_mat_el ) * inv_sqrt_wt(i)       
       DO i = 2, n
          DO j = 1, i
             tr(i,j) = inv_sqrt_wt(i) * tr(i,j) * inv_sqrt_wt(j)
             tr(j,i) = tr(i,j)
          END DO
       END DO
  ELSE
!
!      This is the general case when there are contributions at both ends of the sector.
!      at the left hand boundary.
!       
       i = 1
       tr(i,i) = inv_sqrt_wt(i) * ( tr(i,i) + left_end_mat_el ) * inv_sqrt_wt(i)       
       DO i = 2, n - 1
          DO j = 1, i
             tr(i,j) = inv_sqrt_wt(i) * tr(i,j) * inv_sqrt_wt(j)
             tr(j,i) = tr(i,j)
          END DO
       END DO
       i = n
       DO j = 1, n - 1
          tr(i,j) = inv_sqrt_wt(i) * tr(i,j) * inv_sqrt_wt(j)
          tr(j,i) = tr(i,j)
       END DO
       tr(i,i) = inv_sqrt_wt(i) * ( tr(i,i) + right_end_mat_el ) * inv_sqrt_wt(i)
  END IF
END SUBROUTINE Re_KE
