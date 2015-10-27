!deck Matrix_Renormalization.f
!***begin prologue     Matrix_Renormalization
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
!***end prologue       Matrix_Renormalization

  SUBROUTINE Matrix_Renormalization(grid,eo)
  IMPLICIT NONE
  TYPE(regional_grid)       :: grid  
  INTEGER                   :: eo
  REAL*8                    :: dum
!
!

  IF ( nreg == 1 ) THEN
!
!      One Region only.  No change needed except normalization.
!
       i = 1
       Call Re_KE (grid%reg_mat(i,eo)%tr,                                     &
                   dum,                                                       &
                   dum,                                                       &
                   grid%reg_pt_wt(i)%inv_sqrt_wtr,                            &
                   npt(i),                                                    &
                   i )
  ELSE
       i = 1
!
!      Region one is a special case.  Only the last right hand element needs to be changed
!      and the elements renormalized.
!
       Call Re_KE (grid%reg_mat(i,eo)%tr,                                     &
                   dum,                                                       &
                   reg_grid(icoord)%reg_mat(i+1,eo)%tr(1,1),                  &
                   reg_grid(icoord)%reg_pt_wt(i)%inv_sqrt_wtr,                &
                   npt(i),                                                    &
                   i )
!
!      Now do the general case.  Both the left and right elemenents are changed.
! 
       DO i = 2 , nreg - 1
          Call Re_KE (reg_grid(icoord)%reg_mat(i,eo)%tr,                      &
                      reg_grid(icoord)%reg_mat(i-1,eo)%tr(npt(i-1),npt(i_1)), &
                      reg_grid(icoord)%reg_mat(i+1,eo)%tr(1,1),               &
                      reg_grid(icoord)%reg_pt_wt(i)%inv_sqrt_wtr,             &
                      npt(i),                                                 &
                      i )
       END DO
       i = nreg
!
!      Special case of last region.  Only the left hand elements needs to be changed.
!
          Call Re_KE (reg_grid(icoord)%reg_mat(i,eo)%tr,                      &
                      reg_grid(icoord)%reg_mat(i-1,eo)%tr(npt(i-1),npt(i_1)), &
                      dum,                                                    &
                      reg_grid(icoord)%reg_pt_wt(i)%inv_sqrt_wtr,             &
                      npt(i),                                                 &
                      i )
  END IF
END SUBROUTINE Matrix_Renormalization
