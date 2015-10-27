!deck Ke_Fedvr.f
!***begin prologue     Ke_Fedvr
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
!***end prologue       Ke_Fedvr

  SUBROUTINE Ke_Fedvr(type)
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  CHARACTER(LEN=*)                     :: type
  REAL*8                               :: dum
  REAL*8, DIMENSION(:), ALLOCATABLE    :: fac
  REAL*8                               :: one = 1.d0
  INTEGER                              :: icoord
  INTEGER                              :: max_val
!
  max_val = 0
  DO i = 1, nreg
     max_val = mAx(max_val,npt(i))
  END DO
  ALLOCATE ( fac( 1 : max_val ) )
  IF ( type == 'eta') THEN
       icoord = 1
       DO i = nreg
          ALLOCATE( reg_grid(icoord)%reg_mat(i,0)%(tr( npt(i), npt(i) ) )
          fac( 1 : npt(i) ) = one - reg_grid(icoord)%reg_pt_wt(i)%qr( 1 : npt(i) )         &
                                            *                                              &
                                   reg_grid(icoord)%reg_pt_wt(i)%qr( 1 : npt(i) )
          CALL Eta_KE_Even ( reg_grid(icoord)%reg_mat(i,0)%tr,                             &
                             reg_grid(icoord)%reg_pt_wt(i)%qr,                             &
                             reg_grid(icoord)%reg_pt_wt(i))%wtr,                           &
                             reg_grid(icoord)%reg_poly(i,0)%pr,                            &
                             reg_grid(icoord)%reg_poly(i,0)%dpr,                           &
                             fac,                                                          &
                             npt(i) )
          IF (m_max > 0 ) THEN
              ALLOCATE( reg_grid(icoord)%reg_mat(i,1)%(tr( npt(i), npt(i) ) )
              CALL Eta_KE_Odd  ( reg_grid(icoord)%reg_mat(i,1)%tr,                         &
                                 reg_grid(icoord)%reg_pt_wt(i)%qr,                         &
                                 reg_grid(icoord)%reg_pt_wt(i))%wtr,                       &
                                 reg_grid(icoord)%reg_poly(i,1)%pr,                        &
                                 reg_grid(icoord)%reg_poly(i,1)%dpr,                       &
                                 fac,                                                      &                                               
                                 npt(i) )                                               
          END IF
       END DO
       Call Matrix_Renormalization(icoord,0)
       IF (m_max > 0 ) THEN       
           Call Matrix_Renormalization(icoord,1)
       END IF
  ELSE IF ( type == 'xi' ) THEN
       icoord = 2
       DO i = nreg
          ALLOCATE( reg_grid(icoord)%reg_mat(i,0)%(tr( npt(i), npt(i) ) )
          fac( 1 : npt(i) ) = reg_grid(icoord)%reg_pt_wt(i)%qr( 1 : npt(i) )               &
                                            *                                              &
                             reg_grid(icoord)%reg_pt_wt(i)%qr( 1 : npt(i) )                &
                                            -                                              &
                                           one
          CALL Xi_KE_Even  ( reg_grid(icoord)%reg_mat(i,0)%tr,                             &
                             reg_grid(icoord)%reg_pt_wt(i)%qr,                             &
                             reg_grid(icoord)%reg_pt_wt(i))%wtr,                           &
                             reg_grid(icoord)%reg_poly(i,0)%pr,                            &
                             reg_grid(icoord)%reg_poly(i,0)%dpr,                           &
                             fac,                                                          &
                             npt(i) )
          IF (m_max > 0 ) THEN
              ALLOCATE( reg_grid(icoord)%reg_mat(i,1)%(tr( npt(i), npt(i) ) )
              CALL Xi_KE_Odd   ( reg_grid(icoord)%reg_mat(i,1)%tr,                         &
                                 reg_grid(icoord)%reg_pt_wt(i)%qr,                         &
                                 reg_grid(icoord)%reg_pt_wt(i))%wtr,                       &
                                 reg_grid(icoord)%reg_poly(i,1)%pr,                        &
                                 reg_grid(icoord)%reg_poly(i,1)%dpr,                       &
                                 fac,                                                      &
                                 npt(i) )
          END IF
       END DO
       Call Matrix_Renormalization(icoord,0)
       IF (m_max > 0 ) THEN       
           Call Matrix_Renormalization(icoord,1)
       END IF
  END IF
  DEALLOCATE ( fac )
!
END SUBROUTINE Ke_Fedvr
