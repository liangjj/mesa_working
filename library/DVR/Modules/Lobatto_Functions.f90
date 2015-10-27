!deck Lobatto_Functions.f
!***begin prologue     Lobatto_Functions
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
!***end prologue       Lobatto_Functions

  SUBROUTINE Lobatto_Functions(typwt,coord)
  USE dvr_global
  IMPLICIT NONE
  CHARACTER(LEN=*)          :: typwt
  CHARACTER(LEN=*)          :: coord
  REAL*8                    :: one = 1.d0
  REAL*8                    :: dum
  INTEGER                   :: i
  INTEGER                   :: j
  INTEGER                   :: icoord
!
!
  Write(iout,1) coord
  IF (coord == 'eta' ) THEN
      icoord = 1
  ELSE IF (coord == 'xi' ) THEN
      icoord = 2
  END IF
  DO  i=1,nreg
      Write(iout,2) i, npt(i)
!
!        calculate the sector functions and their derivatives for the m even grid.
!
      ALLOCATE( reg_grid(icoord)%reg_pt_wt(i)%qr(npt(i)),                                 &
                reg_grid(icoord)%reg_pt_wt(i)%wtr(npt(i)),                                &
                reg_grid(icoord)%reg_poly(i,0)%pr(npt(i),npt(i)),                         &
                reg_grid(icoord)%reg_poly(i,0)%dpr(npt(i),npt(i)),                        &
                reg_grid(icoord)%reg_poly(i,0)%ddpr(npt(i),npt(i)),                       &
                reg_grid(icoord)%reg_pt_wt(i)%inv_sqrt_wtr(npt(i)) )              
      CALL drvply( reg_grid(icoord)%reg_pt_wt(i)%qr,                                      &
                   reg_grid(icoord)%reg_pt_wt(i)%wtr,                                     &
                   reg_grid(icoord)%reg_poly(i,0)%pr,                                     &
                   reg_grid(icoord)%reg_poly(i,0)%dpr,                                    &
                   reg_grid(icoord)%reg_poly(i,0)%ddpr,                                   &
                   edge(i),                                                               &
                   typwt,                                                                 &
                   npt(i),                                                                &
                   i)
  END DO
!
!                       To compute most of what is required, it is not necessary
!                       to construct anything else. Since the bridge functions 
!                       span two elements, one can define them at the grid points
!                       but their derivatives are discontinuous across the
!                       sector boundaries. The matrix elements can be constructed 
!                       entirely from re-normalized sector quantities. 
!
  IF ( nreg == 1) THEN
!
!                       Only one region.  No endpoint corrections required.
!
       i = 1
       Call ReGrid ( reg_grid(icoord)%reg_pt_wt(i)%wtr,                                   &
                     reg_grid(icoord)%reg_pt_wt(i)%inv_sqrt_wtr,                          &
                     dum,                                                                 &
                     dum,                                                                 &
                     npt(i),                                                              &
                     i)               
  ELSE
!
!                       First region.  Correction  at right endpoint needed from first function
!                       in refion 2.
       i = 1
       Call ReGrid ( reg_grid(icoord)%reg_pt_wt(i)wtr,                                    &
                     reg_grid(icoord)%reg_pt_wt(i)%inv_sqrt_wtr,                          &
                     dum,                                                                 &
                     reg_grid(icoord)%reg_pt(i+1)%wtr(1),                                 &
                     npt(i),                                                              &
                     i)        
!       
       DO i = 2, nreg - 1
!
!                       General case.  Put result from the previous region into the
!                       the left region and correct the right endpoint.
!
          Call ReGrid ( reg_grid(icoord)%reg_pt_wt(i)%wtr,                                 &
                        reg_grid(icoord)%reg_pt_wt(i)%inv_sqrt_wtr,                        &
                        reg_grid(icoord)%reg_pt_wt(i-1)%wtr(npt(i-1)),                     &
                        reg_grid(icoord)%reg_pt_wt(i+1)%wtr(1),                            &
                        npt(i),                                                            &
                        i)               
       END DO
!
!                       Last region.  Correct the left end point.
!
       i = nreg
          Call ReGrid ( reg_grid(icoord)%reg_pt_wt(i)%wtr,                                 &
                        reg_grid(icoord)%reg_pt_wt(i)%inv_sqrt_wtr,                        &
                        reg_grid(icoord)%reg_pt_wt(i-1)%wtr(npt(i-1)),                     &
                        dum,                                                               &
                        npt(i),                                                            &
                        i)               
  END IF
  DO i = 1, nreg 
     Call RePoly (reg_grid(icoord)%reg_poly(i,0)%pr,                                       &
                  reg_grid(icoord)%reg_poly(i,0)%dpr,                                      &
                  reg_grid(icoord)%reg_poly(i,0)%ddpr,                                     &
                  reg_grid(icoord)%reg_pt_wt(i)%inv_sqrt_wtr,                              &
                  npt(i),                                                                  &
                  i)             
  END DO
!     
  IF ( m_max >= 0 ) THEN
!
!         Do the same for the m odd grid.
!
       DO i = 1,nreg
          ALLOCATE( reg_grid(icoord)%reg_poly(i,1)%pr(npt(i),npt(i)),                      &
                    reg_grid(icoord)%reg_poly(i,1)%dpr(npt(i),npt(i)),                     &
                    reg_grid(icoord)%reg_poly(i,1)%ddpr(npt(i),npt(i)) )
          Call renorm( reg_grid(icoord)%reg_pt_wt(i)%qr,                                   &
                       reg_grid(icoord)%reg_poly(i,0)%pr,                                  & 
                       reg_grid(icoord)%reg_poly(i,0)%dpr,                                 & 
                       reg_grid(icoord)%reg_poly(i,0)%ddpr,                                & 
                       reg_grid(icoord)%reg_poly(i,1)%pr,                                  & 
                       reg_grid(icoord)%reg_poly(i,1)%dpr,                                 & 
                       reg_grid(icoord)%reg_poly(i,1)%ddpr,                                & 
                       npt(i),                                                             &
                       coord )                               
       END DO
!
! 
                      Now make the proper normalized quantities
!
      DO i = 1, nreg
         Call RePoly ( reg_grid(icoord)%reg_poly(i,1)%pr,                                  &
                       reg_grid(icoord)%reg_poly(i,1)%dpr,                                 & 
                       reg_grid(iccord)%reg_poly(i,1)%ddpr,                                &
                       reg_grid(icoord)%pt_wt(i)%inv_sqrt_wtr,                             &
                       npt(i),                                                             &
                       i)             
      END DO
  END IF
1 FORMAT(/,10x,'Calculate the Regional Basis Functions for Coordinate = ',a4)
2 FORMAT(/,10x,'Region = ',i3,2x,'Number of Lobatto_Functions = ', i3)
END SUBROUTINE Lobatto_Functions
