!deck H_0_Fedvr.f
!***begin prologue     H_0_Fedvr
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate the final unperturbed Hamiltonian for each m value.
!***                   
!***                   
!***                   
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       H_0_Fedvr

  SUBROUTINE H_0_Fedvr(type)
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  CHARACTER(LEN=*)                    :: type
  REAL*8, DIMENSION(:), ALLOCATABLE   :: fac
  INTEGER                             :: icoord
  INTEGER                             :: max_val
  INTEGER                             :: i
  INTEGER                             :: j
!
!
  max_val = 0
  DO i = 1, nreg
     max_val = max ( max_val, npt(i) )
  END DO
  ALLOCATE( fac( 1 : max_val ) )
  IF ( type == 'eta') THEN
       icoord = 1
       DO i = 1, nreg
          fac ( 1 : npt(i) ) = 1. d0 / ( one - reg_grid(icoord)%reg_pt_wt(i)%qr ( 1 : npt(i) )       &
                                          *                                                          &
                                               reg_grid(icoord)%reg_pt_wt(i)%qr ( 1 : npt(i) ) )
          ALLOCATE( reg_grid(icoord)%reg_mat(i,0)%(ham( 1 :  npt(i), 1 : npt(i) ) )
          reg_grid(icoord)%reg_mat(i,0)%ham ( :,: ) = reg_grid(icoord)%reg_mat(i,0)%tr ( :,: ) 
          DO j = 1 , npt(i)
             reg_grid(icoord)%reg_mat(i,0)%ham ( j,j )                                               &
                                            =                                                        &
             reg_grid(icoord)%reg_mat(i,0)%ham ( j,j )                                               &
                                            +                                                        &
                                            reg_grid(icoord)%reg_mat(i,0)%vr ( j )         
          END DO                 
          DO m = 2, m_max,2
             ALLOCATE( reg_grid(icoord)%reg_mat(i,m)%(ham( 1 :  npt(i), 1 : npt(i) ) )
             reg_grid(icoord)%reg_mat(i,m)%ham ( : , : )                         )                   &
                                          =                                                          &
             reg_grid(icoord)%reg_mat(i,0)%ham ( : , : ) 
             DO j = 1, npt(i)
                reg_grid(icoord)%reg_mat(i,m)%ham ( j ,j ) 
                                            =                                                        &
                reg_grid(icoord)%reg_mat(i,m)%ham ( j , j )                                          &
                                            -                                                        &
                                              m * m * fac ( j )                 
             END DO
          END DO
          IF ( m_max > 0 ) THEN
               ALLOCATE( reg_grid(icoord)%reg_mat(i,1)%(ham( 1 :  npt(i), 1 : npt(i) ) )
               reg_grid(icoord)%reg_mat(i,1)%ham ( :,: ) = reg_grid(icoord)%reg_mat(i,1)%tr ( :,: ) 
               DO j = 1 , npt(i)
                  reg_grid(icoord)%reg_mat(i,1)%ham ( j,j )                                          &
                                                 =                                                   &
                  reg_grid(icoord)%reg_mat(i,1)%ham ( j,j )                                          &
                                                 +                                                   &
                                                 reg_grid(icoord)%reg_mat(i,0)%vr ( j )         
               END DO                 
               DO m = 3, m_max,2
                 ALLOCATE( reg_grid(icoord)%reg_mat(i,m)%(ham( 1 :  npt(i), 1 : npt(i) ) )
                 reg_grid(icoord)%reg_mat(i,m)%ham ( 1 : npt(i), 1 : npt(i) )                        &
                                             =                                                       &
                 reg_grid(icoord)%reg_mat(i,0)%ham ( 1 : npt(i), 1 : npt(i) ) 
                 DO j = 1, npt(i)
                    reg_grid(icoord)%reg_mat(i,m)%ham ( j , j )                                      &
                                                =                                                    &
                    reg_grid(icoord)%reg_mat(i,m)%ham ( j , j )                                      &
                                                      -                                              &
                                                  m * m * fac ( j )                 
                 END DO
               END DO
          END IF
          DEALLOCATE( fac )
       END DO
!
  ELSE IF ( type == 'xi' ) THEN
       icoord = 2
       DO i = nreg
          ALLOCATE( reg_grid(icoord)%reg_mat(i,0)%(tr( npt(i), npt(i) ) )
          CALL Xi_KE_Even  ( reg_grid(icoord)%reg_mat(i,0)%tr,                              &
                             reg_grid(icoord)%reg_pt_wt(i)%qr,                             &
                             reg_grid(icoord)%reg_pt_wt(i))%wtr,                           &
                             reg_grid(icoord)%reg_poly(i,0)%pr,                            &
                             reg_grid(icoord)%reg_poly(i,0)%dpr,                           &
                             npt(i),                                                       &
                             i)
       END DO
       Call Matrix_Renormalization(icoord,0)
       IF (m_max > 0 ) THEN
           DO i = nreg
              ALLOCATE( reg_grid(icoord)%reg_mat(i,1)%(tr( npt(i), npt(i) ) )
              CALL Xi_KE_Odd   ( reg_grid(icoord)%reg_mat(i,0)%tr,                          &
                                 reg_grid(icoord)%reg_pt_wt(i)%qr,                         &
                                 reg_grid(icoord)%reg_pt_wt(i))%wtr,                       &
                                 reg_grid(icoord)%reg_poly(i,1)%pr,                        &
                                 reg_grid(icoord)%reg_poly(i,1)%dpr,                       &
                                 npt(i),                                                   &
                                 i)
           END DO
          Call Matrix_Renormalization(icoord,1)
       END IF
  END IF
!
END SUBROUTINE H_0_Fedvr
