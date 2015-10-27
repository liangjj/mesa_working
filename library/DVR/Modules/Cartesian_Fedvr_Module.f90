\!***********************************************************************
! Cartesian_Fedvr_Module
!**begin prologue     Cartesian_Fedvr_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***description       
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                  
!***                
!***       
!***               
!***
!***               
!***               
!                  
!                  
!                  
!                  
!                  
!***               
!***               
!***               
!***               
!***
!***              
!***!***          
!***              
!***              
!***              
!***references
!***modules needed    See USE statements below
!***comments      
!***              
!***              
!***              
!***              
!***end prologue      Cartesian_Fedvr_Module
!***********************************************************************
!***********************************************************************
                           MODULE Cartesian_Fedvr_Module
                           USE Spheroidal_DVR_Global
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
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

  SUBROUTINE Lobatto_Functions(grid,typwt,keyword)
  IMPLICIT NONE
  TYPE(regional_grid)                  :: grid
  CHARACTER(LEN=*)                     :: typwt
  CHARACTER(LEN=*)                     :: keyword
  REAL*8                               :: one = 1.d0
  REAL*8                               :: dum
  INTEGER                              :: i
!
!
  ALLOCATE ( grid%reg_pt_wt(1:nreg), grid%reg_poly(1:nreg) )
  Call Polynomials(grid,grid%reg_poly,typwt)
1 FORMAT(/,10x,'Calculate the Regional Basis Functions for Coordinate = ',a4)
2 FORMAT(/,10x,'Region = ',i3,2x,'Number of Polynomials = ', i3)
END SUBROUTINE Lobatto_Functions
!***********************************************************************
!***********************************************************************
!deck Coordinate_Factors.f
!***begin prologue     Coordinate_Factors
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Coordinate_Factors(grid,keyword)
  IMPLICIT NONE
  TYPE(regional_grid)                  :: grid
  CHARACTER (LEN=*)                    :: keyword
  INTEGER                              :: i
!
  IF ( keyword == 'eta') THEN
       DO i = 1, nreg
          ALLOCATE ( grid%reg_pt_wt(i)%qr_fac( 1:npt(i) ),           &
                     grid%reg_pt_wt(i)%inv_qr_fac( 1:npt(i) ),       &
                     grid%reg_pt_wt(i)%inv_sqrt_qr_fac( 1:npt(i) ) )
          grid%reg_pt_wt(i)%qr_fac(:) = one - grid%reg_pt_wt(i)%qr(:) * grid%reg_pt_wt(i)%qr(:)
          grid%reg_pt_wt(i)%inv_qr_fac(:) = ( one / grid%reg_pt_wt(i)%qr_fac(:) )
          grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = Sqrt ( grid%reg_pt_wt(i)%inv_qr_fac(:) )
       END DO
  ELSE IF ( keyword == 'xi') THEN
       DO i = 1, nreg
          ALLOCATE ( grid%reg_pt_wt(i)%qr_fac( 1:npt(i) ),           &
                     grid%reg_pt_wt(i)%inv_qr_fac( 1:npt(i) ),       &
                     grid%reg_pt_wt(i)%inv_sqrt_qr_fac( 1:npt(i) ) )
          grid%reg_pt_wt(i)%qr_fac(:) = grid%reg_pt_wt(i)%qr(:) * grid%reg_pt_wt(i)%qr(:) - one
          grid%reg_pt_wt(i)%inv_qr_fac(:) = ( one / grid%reg_pt_wt(i)%qr_fac(:) )
          grid%reg_pt_wt(i)%inv_sqrt_qr_fac(:) = Sqrt ( grid%reg_pt_wt(i)%inv_qr_fac(:) )
       END DO
  END IF
!
!
END SUBROUTINE Coordinate_Factors
!***********************************************************************
!***********************************************************************
!deck FE_DVR_Matrices.f
!***begin prologue     FE_DVR_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            
!***                   
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE FE_DVR_Matrices(grid)
  IMPLICIT NONE
  TYPE(regional_grid)                  :: grid
!
!
  ALLOCATE(reg_mat_e(1:nreg))
  Call KE_FEDVR_Even(grid,grid%reg_mat_e)
  IF ( m_max > 0 ) THEN
       ALLOCATE(reg_mat_o(1:nreg))
       Call KE_FEDVR_Even(grid,grid%reg_mat_o)
  END IF
!
!
END SUBROUTINE FE_DVR_Matrices
!***********************************************************************
!***********************************************************************
!deck KE_FEDVR_Even.f
!***begin prologue     KE_FEDVR_Even
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
!***end prologue       Ke_Fedvr_Even

  SUBROUTINE KE_FEDVR_Even(grid,reg_mat_e)

  IMPLICIT NONE
  TYPE(regional_grid)                  :: grid  
  TYPE(even_matrices), DIMENSION(:)    :: reg_mat_e  
  INTEGER                              :: ireg
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: k
!
!
  DO ireg = 1, nreg
     ALLOCATE( grid%reg_mat(ireg,0)%tr( npt(ireg), npt(ireg) ) )
     grid%reg_mat(ireg,0)%tr(:,:) = zero
     DO i = 1, npt(ireg)
        DO j = 1, i
           DO k = 1, npt(ireg)
              grid%reg_mat(ireg,0)%tr(i,j)                               &
                                =                                        &
              grid%reg_mat(ireg,0)%tr(i,j)                               &
                                -                                        &
              grid_pt_wt(ireg)%ptr_fac(k)                                &
                                *                                        &
              grid%reg_pt_wt(ireg)%wtr(k)                                & 
                                *                                        &
              grid%reg_poly(ireg,0)%dpr(k,i)                             & 
                                *                                        &
              grid%reg_poly(ireg,0)%dpr(k,j) 
        END DO
        grid%reg_mat(ireg,0)%tr(j,i) = grid%reg_mat(ireg,0)%tr(i,j)  
     END DO
  END DO
  Call Matrix_Renormalization(grid,0)
!
END SUBROUTINE Ke_Fedvr_Even
!***********************************************************************
!***********************************************************************
!deck Ke_Fedvr_Odd.f
!***begin prologue     Ke_Fedvr_Odd
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
!***end prologue       Ke_Fedvr_Odd

  SUBROUTINE KE_FEDVR_Even(grid,reg_mat_o)
  IMPLICIT NONE
  TYPE(regional_grid)                  :: grid  
  TYPE(even_matrices), DIMENSION(:)    :: reg_mat_o  
  INTEGER                              :: ireg
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: k
!
!
  DO ireg = 1, nreg
     ALLOCATE( grid%reg_mat(ireg,1)%tr( npt(ireg), npt(ireg) ) )
     grid%reg_mat(ireg,1)%tr(:,:) = zero
     DO i = 1, npt(ireg)
        DO j = 1, i
           DO k = 1, npt(ireg)
              grid%reg_mat(ireg,1)%tr(i,j)                               &
                                =                                        &
              grid%reg_mat(ireg,1)%tr(i,j)                               &
                                -                                        &
              grid_reg_pt_wt(ireg)%ptr_fac(k)                            &
                                *                                        &
              grid_reg_pt_wt(ireg)%ptr_fac(k)                            &
                                *
              grid%reg_pt_wt(ireg)%wtr(k)                                & 
                                *                                        &
              grid%reg_poly(ireg,1)%dpr(k,i)                             & 
                                *                                        &
              grid%reg_poly(ireg,1)%dpr(k,j) 
           END DO 
           grid%reg_mat(ireg,1)%tr(i,j)                                  &
                             =                                           &
           grid%reg_mat(ireg,1)%tr(i,j)                                  &
                             +                                           &
           grid%reg_pt_wt(ireg)%qr(i)                                    &
                             *                                           &
           grid_reg_pt_wt(ireg)%ptr_fac(i)                               &
                             *                                           &
           grid%reg_poly(ireg,1)%pr(i,i)                                 &
                             *                                           &
           grid%reg_poly(ireg,1)%dpr(i,j)                                &
                             +                                           &
           grid%reg_pt_wt(ireg)%qr(j)                                    &
                             *                                           &
           grid_reg_pt_wt(ireg)%ptr_fac(j)                               &
                             *                                           &
           grid%reg_poly(ireg,1)%pr(j,j)                                 &
                             *                                           &
           grid%reg_poly(ireg,1)%dpr(j,i)
!
           grid%reg_mat(ireg,1)%tr(j,i) =  grid%reg_mat(ireg,1)%tr(i,j)
     END DO
     grid%reg_mat(ireg,1)%tr(i,i) =  grid%reg_mat(ireg,1)%tr(i,i)        &
                                  -                                      &
                         grid%reg_pt_wt(ireg)%qr(i)                      &
                                  *                                      &
                         grid%reg_pt_wt(ireg)%qr(i)                      &
                                  *                                      &
                         grid%reg_poly(ireg,1)%pr(i)                     & 
                                  *                                      &
                         grid%reg_poly(ireg,1)%pr(i)                     & 
  END DO
  Call Matrix_Renormalization (grid,1)
!
END SUBROUTINE Ke_Fedvr_Odd
!***********************************************************************
!***********************************************************************
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
          DO i = 1 , npt(i)
             reg_grid(icoord)%reg_mat(i,0)%tr ( i,i ) 
                                            =                                              &
             reg_grid(icoord)%reg_mat(i,0)%ham ( i,i )                                     &
                                            +                                              &
                                            reg_grid(icoord)%reg_mat(i,0)%vr ( i )         
          END DO                 
          DO m = 2, m_max,2
             ALLOCATE( reg_grid(icoord)%reg_mat(i,m)%(ham( 1 :  npt(i), 1 : npt(i) ) )
             reg_grid(icoord)%reg_mat(i,m)%ham ( 1 : npt(i), 1 : npt(i) )                   &
                                          =                                                 &
             reg_grid(icoord)%reg_mat(i,0)%ham ( 1 : npt(i), 1 : npt(i) ) 
             DO i = 1, npt(i)
                reg_grid(icoord)%reg_mat(i,m)%ham ( i,i ) 
                                            =                                              &
                reg_grid(icoord)%reg_mat(i,m)%ham ( i,i )                                  &
                                            -                                              &
                                              m * m * fac ( i )                 
             END DO
          END DO
          DEALLOCATE( fac )
       END DO
                  
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
              CALL Eta_KE_Odd  ( reg_grid(icoord)%reg_mat(i,0)%tr,                          &
                                 reg_grid(icoord)%reg_pt_wt(i)%qr,                         &
                                 reg_grid(icoord)%reg_pt_wt(i))%wtr,                       &
                                 reg_grid(icoord)%reg_poly(i,1)%pr,                        &
                                 reg_grid(icoord)%reg_poly(i,1)%dpr,                       &
                                 npt(i),                                                   &
                                 i)
           END DO
       END IF
       Call Matrix_Renormalization(icoord,1)
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
!***********************************************************************
!***********************************************************************

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

  SUBROUTINE Pe_Fedvr(type)
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  CHARACTER(LEN=*)          :: type
  REAL*8                    :: dum
  INTEGER                   :: icoord
!
!
  IF ( type == 'eta') THEN
       icoord = 1
       DO i = nreg
          ALLOCATE( reg_grid(icoord)%reg_mat(i,0)%vr( npt(i) ) )
          reg_grid(icoord)%reg_mat(i,0)%vr( : ) = R_ab * ( z_b - z_a )                     &
                                                       *                                   &
                                                   reg_grid(icoord)%reg_pt_wt(i)%qr( : )
       END DO
!
  ELSE IF ( type == 'xi' ) THEN
       icoord = 2
       DO i = nreg
          ALLOCATE( reg_grid(icoord)%reg_mat(i,0)%vr( npt(i) ) )
          reg_grid(icoord)%reg_mat(i,0)%vr( : ) = R_ab * ( z_a + z_b )                     &
                                                       *                                   &
                                                   reg_grid(icoord)%reg_pt_wt(i)%qr( : )
       END DO
  END IF
!
END SUBROUTINE Pe_Fedvr
!***********************************************************************
!***********************************************************************
!deck Renorm.f
!***begin prologue     Renorm
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
!***end prologue       Renorm

  SUBROUTINE Renorm(q,p_e,dp_e,ddp_e,p_o,dp_o,ddp_o,n,coord)
  IMPLICIT NONE
  REAL*8, DIMENSION (:)     :: q
  REAL*8, DIMENSION (:,:)   :: p_e
  REAL*8, DIMENSION (:,:)   :: dp_e
  REAL*8, DIMENSION (:,:)   :: ddp_e
  REAL*8, DIMENSION (:,:)   :: p_o
  REAL*8, DIMENSION (:,:)   :: dp_o
  REAL*8, DIMENSION (:,:)   :: ddp_o
  CHARACTER(LEN=*)          :: coord
  REAL*8                    :: one = 1.d0
  REAL*8                    :: fac
  INTEGER                   :: i
!
!
  IF (coord == 'eta' ) THEN
      DO i = 1, n
         fac = Sqrt ( one / ( one - q(i) * q(i) ) )
         p_o(:,i)  = fac * p_e(:,j) 
         dp_o(:,i) = fac * dp_e(:,i) 
         ddp_o(:,i) = fac * ddp_e(:,i) 
      END DO
  ELSE IF (coord == 'xi' ) THEN
      DO i= 1, n
         fac = Sqrt ( one / ( q(i) * q(i) - one ) )
         p_o(:,i) = fac * p_e(:,i) 
         dp_o(:,i) = fac * dp_e(:,i) 
         ddpr_o(:,i) = fac * ddp_e(:,i) 
      END DO
  END IF
END SUBROUTINE Renorm
!***********************************************************************
!***********************************************************************
!deck Eta_Coordinate_Factors                                                  
!***begin prologue     Eta_Coordinate_Factors
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author
!***source                                                                                                           
!***purpose                                                                                                                             
!***                                                                                                                                   
!***                                                                                                                      
!***description                                                                                                                    
!***                                                                                                                                    
!***                                                                                                                                     
!***                                                                                                                                    
!***references                                                                                                                          
!***routines called                                                                                                                      
!***end prologue       Eta_Coordinate_Factors                                                                                            
  Subroutine Eta_Coordinate_Factors (grid)
  IMPLICIT NONE
  TYPE(regional_grid)                  :: grid
  REAL*8                               :: one = 1.d0
  INTEGER                              :: i
  DO i = 1, nreg
     ALLOCATE ( grid%reg_pt_wt(i)%ptr_fac( 1:npt(i) ),           &
                grid%reg_pt_wt(i)%inv_ptr_fac( 1:npt(i) ),       &
                grid%reg_pt_wt(i)%inv_sqrt_ptr( 1:npt(i) ) )
     grid%reg_pt_wt(i)%ptr_fac(:) = one - grid%reg_pt_wt(i)%ptr(:) * grid%reg_pt_wt(i)%ptr(:)
     grid%reg_pt_wt(i)%inv_ptr_fac(:) = ( one / grid%reg_pt_wt(i)%ptr_fac(:) )
     grid%reg_pt_wt(i)%inv_sqrt_ptr(:) = Sqrt ( grid%reg_pt_wt(i)%inv_ptr_fac(:) )
  END DO
  DEALLOCATE ( tmp )
  END Subroutine Eta_Coordinate_Factors
!***********************************************************************
!***********************************************************************
!deck Xi_Coordinate_Factors
!***begin prologue     Xi_Coordinate_Factors
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author
!***source
!***purpose                                                                                                                             
!***
!***description
!***                                                                                                            
!***                                                                                                                 
!***references                                                                                                           
!***routines called                                                                                                          
!***end prologue       Xi_Coordinate_Factors                                                                                     
  Subroutine Xi_Coordinate_Factors (grid)
  IMPLICIT NONE
  TYPE(regional_grid)                  :: grid
  REAL*8                               :: one = 1.d0
  INTEGER                              :: i
  DO i = 1, nreg
     ALLOCATE ( grid%reg_pt_wt(i)%ptr_fac( 1:npt(i) ),           &
                grid%reg_pt_wt(i)%inv_ptr_fac( 1:npt(i) ),       &
                grid%reg_pt_wt(i)%inv_sqrt_ptr( 1:npt(i) ) )
     grid%reg_pt_wt(i)%ptr_fac(:) = grid%reg_pt_wt(i)%ptr(:) * grid%reg_pt_wt(i)%ptr(:) - one
     grid%reg_pt_wt(i)%inv_ptr_fac(:) = ( one / grid%reg_pt_wt(i)%ptr_fac(:) )
     grid%reg_pt_wt(i)%inv_sqrt_ptr(:) = Sqrt ( grid%reg_pt_wt(i)%inv_ptr_fac(:) )
  END DO
  DEALLOCATE ( tmp )
  END Subroutine Xi_Coordinate_Factors
!***********************************************************************
!***********************************************************************
           END MODULE Spheroidal_Fedvr_Module
!***********************************************************************
!***********************************************************************
