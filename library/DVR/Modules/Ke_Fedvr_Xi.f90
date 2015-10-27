!deck Ke_Fedvr_Xi.f
!***begin prologue     Ke_Fedvr_Xi
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
!***end prologue       Ke_Fedvr_Xi

  SUBROUTINE Ke_Fedvr_Xi
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  REAL*8                    :: dum
!
!
  DO  i=1,nreg
!
!
!        Calculate the kinetic energy matrix.
!        this is done over un-normalized functions.  Renormalization
!        is necessary when joining different sectors.
!
      ALLOCATE( reg_eta_mat(i,2,1)%tr( npt(i), npt(i) ) )
!
      CALL Xi_KE ( reg_eta_mat(i,1,1)%tr, reg_grid(i)%qr, reg_grid(i)%wtr,                               &
                    reg_grid(i)%pr_e, reg_grid(i)%dpr_e, npt(i), i)
  END DO
  IF ( nreg == 1 ) THEN
!
!      One Region only.  No change needed except normalization.
!
       i = 1
       Call Re_KE (reg_eta_mat(i,2,1)%tr, dum, dum, reg_grid(i)%inv_sqrt_wtr, npt(i), i )
  ELSE
       i = 1
!
!      Region one is a special case.  Only the last right hand element needs to be changed
!      and the elements renormalized.
!
       Call Re_KE (reg_eta_mat(i,2,1)%tr, dum, reg_eta_mat(i+1,2,1)%t(1,1),                               &
                   reg_grid(i)%inv_sqrt_wtr, npt(i), i )
!
!      Now do the general case.  Both the left and right elemenents are changed.
! 
       DO i = 2 , nreg - 1
          Call Re_KE (reg_eta_mat(i,2,1)%tr, reg_eta_mat(i-1,2,1)%tr(npt(i-1),npt(i-1)),                  &
                      reg_eta_mat(i+1)%tr(1,1), reg_grid(i)%inv_sqrt_wtr, npt(i), i )
       END DO
       i = nreg
!
!      Special case of last region.  Only the left hand elements needs to be changed.
!
       Call Re_KE (reg_eta_mat(i,2,1)%tr, reg_eta_mat(i-1)%tr(npt(i-1),npt(i-1), dum,                     &
                   reg_grid(i)%inv_sqrt_wtr, npt(i), i )
  END IF
END SUBROUTINE Ke_Fedvr_Xi
