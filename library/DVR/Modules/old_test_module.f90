!***********************************************************************
! Spheroidal_Fedvr_Module
!**begin prologue     Spheroidal_Fedvr_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains all of the major subroutines to propagate
!***                  a wavefunction in time using the Lanczos
!***                  algorithm.  The routine works for stnadard and generalized 
!***                  symmetric problems.
!***                  Explicit interfaces are used to allow a transparent use of 
!***                  generic subroutines which work for both real and complex vectors.  
!***                  This feature permits a single code to be used for both real and
!***                  imaginary time propagation.
!***description       Given a starting vector, a number of iterations
!***                  are performed until the time propagated solution
!***                  satisfies a fixed accuracy criterion.  The criterion
!***                  used depends on whether one is propagating in real or
!***                  imaginary time.  In imaginary time, the number of 
!***                  iterations at a given step is controlled by the 
!***                  convergence of the eigenvalue of interest.  The
!***                  converged vector provides the starting vector for the
!***                  next time step.  For real time, the number of iterations
!***                  at a given step depends on the RMS deviation of the 
!***                  wavefunction.  When that deviation is smaller than some
!***                  preset value, the iterations are terminated.  The 
!***                  converged wavefunction is used as the starting wavefunction
!***                  at the next time step.
!***       
!***                                  Details
!***
!***                  Perform the Lanczos recursion at a given timestep.
!***                  b_{i+1} S |V_{i+1}> = [ H - a_{i} S ] |V_{i}> 
!                                               -   b_{i} S |V_{i-1}>
!                                        T
!                                   S = U U
!                                    -1     -1  -T
!                                   S  = ( U   U   )
!***                  We rewite this as,
!***                                         -T  -1
!***                  b_{i+1} |X_{i+1}> = [ U H U - a_{i} ] |X_{i}> 
!***                                               -   b_{i} |X_{i-1}>
!***                          
!***                                |X_{i}> = U |V_{i}>
!***
!***                  The Lanczos recursion for a generalized eigenvalue problem requires
!***                  the factoization of the S matrix and the solution of two triangular
!***                  linear systems at each iteration ans well as a matrix multiply.
!***                  The multiply and linear system solves are all comparable in flop count.
!***references
!***modules needed    See USE statements below
!***comments          In this portable version I have disabled all unnecessary
!***                  writing to files.  The original Fortran is commented out.
!***                  In addition, there is no option to compute the autocorrelation
!***                  function as this would require reading and manipulating the
!***                  initial state wavefunction from a file.
!***end prologue      Generalized_Lanczos_Module
!***********************************************************************
!***********************************************************************
                           MODULE Spheroidal_Fedvr_Module
                           USE Spheroidal_DVR_Global
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
                            INTERFACE Polynomials                                     
                       MODULE PROCEDURE Even_Polynomials,                   &
                                        Odd_Polynomials  
                            END INTERFACE Polynomials
!
                            INTERFACE Kinetic_Energy                                     
                       MODULE PROCEDURE Even_Kinetic_Energy,                &
                                        Odd_Kinetic_Energy  
                            END INTERFACE Kinetic_Energy
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
  ALLOCATE ( grid%reg_pt_wt(1:nreg), grid%reg_poly_e(1:nreg) )
  Call Polynomials(grid,grid%reg_poly_e,typwt)
  Call Coordinate_Factors(grid,keyword)
!
!
  IF (m_max > 0 ) THEN
      ALLOCATE ( grid%reg_poly_o(1:nreg) )
      Call Polynomials(grid,grid%reg_poly_o,typwt)
  END IF
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
  ALLOCATE(grid%reg_mat_e(1:nreg))
  Call KE_FEDVR_Even(grid,grid%reg_mat_e)
  IF ( m_max > 0 ) THEN
       ALLOCATE(grid%reg_mat_o(1:nreg))
       Call KE_FEDVR_Even(grid,grid%reg_mat_o)
  END IF
!
!
END SUBROUTINE FE_DVR_Matrices
!***********************************************************************
!***********************************************************************
!deck Even_Polynomials.f
!***begin prologue     Even_Polynomials
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
!***end prologue       Even_Polynomials

  SUBROUTINE Even_Polynomials(grid,reg_poly_e,typwt)
  IMPLICIT NONE
  TYPE(regional_grid)                 :: grid
  TYPE(even_functions),  DIMENSION(:) :: reg_poly_e
  CHARACTER(LEN=*)                    :: typwt
  REAL*8                              :: one = 1.d0
  REAL*8                              :: dum
  INTEGER                             :: i
!
!
  DO  i=1,nreg
      Write(iout,2) i, npt(i)
!
!        calculate the sector functions and their derivatives for the m even grid.
!
      ALLOCATE( grid%reg_pt_wt(i)%qr(npt(i)),                                 &
                grid%reg_pt_wt(i)%wtr(npt(i)),                                &
                grid%reg_poly_e(i)%pr(npt(i),npt(i)),                         &
                grid%reg_poly_e(i)%dpr(npt(i),npt(i)),                        &
                grid%reg_poly_e(i)%ddpr(npt(i),npt(i)),                       &
                grid%reg_pt_wt(i)%inv_sqrt_wtr(npt(i)) )              
      CALL drvply( grid%reg_pt_wt(i)%qr,                                      &
                   grid%reg_pt_wt(i)%wtr,                                     &
                   grid%reg_poly_e(i)%pr,                                     &
                   grid%reg_poly_e(i)%dpr,                                    &
                   grid%reg_poly_e(i)%ddpr,                                   &
                   edge(i),                                                   &
                   typwt,                                                         &
                   npt(i),                                                        &
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
       Call ReGrid ( grid%reg_pt_wt(i)%wtr,                                   &
                     grid%reg_pt_wt(i)%inv_sqrt_wtr,                          &
                     dum,                                                     &
                     dum,                                                     &
                     npt(i),                                                  &
                     i)               
  ELSE
!
!                       First region.  Correction  at right endpoint needed from first function
!                       in refion 2.
       i = 1
       Call ReGrid ( grid%reg_pt_wt(i)%wtr,                                   &
                     grid%reg_pt_wt(i)%inv_sqrt_wtr,                          &
                     dum,                                                     &
                     grid%reg_pt_wt(i+1)%wtr(1),                              &
                     npt(i),                                                  &
                     i)        
!       
       DO i = 2, nreg - 1
!
!                       General case.  Put result from the previous region into the
!                       the left region and correct the right endpoint.
!
          Call ReGrid ( grid%reg_pt_wt(i)%wtr,                                 &
                        grid%reg_pt_wt(i)%inv_sqrt_wtr,                        &
                        grid%reg_pt_wt(i-1)%wtr(npt(i-1)),                     &
                        grid%reg_pt_wt(i+1)%wtr(1),                            &
                        npt(i),                                                &
                        i)               
       END DO
!
!                       Last region.  Correct the left end point.
!
       i = nreg
          Call ReGrid ( grid%reg_pt_wt(i)%wtr,                                 &
                        grid%reg_pt_wt(i)%inv_sqrt_wtr,                        &
                        grid%reg_pt_wt(i-1)%wtr(npt(i-1)),                     &
                        dum,                                                   &
                        npt(i),                                                &
                        i)               
  END IF
!
!                        Normalize the functions using the weights.
!
  DO i = 1, nreg 
     Call Re_Poly(                                                             &
                  grid%reg_poly_e(i)%pr,                                       &
                  grid%reg_poly_e(i)%dpr,                                      &
                  grid%reg_poly_e(i)%ddpr,                                     &
                  grid%reg_pt_wt(i)%inv_sqrt_wtr,                              &
                  npt(i) )
  END DO
1 FORMAT(/,10x,'Calculate the Regional Basis Functions for Coordinate = ',a4)
2 FORMAT(/,10x,'Region = ',i3,2x,'Number of Polynomials = ', i3)
END SUBROUTINE Even_Polynomials
!***********************************************************************
!***********************************************************************
!deck Odd_Polynomials.f
!***begin prologue     Odd_Polynomials
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
!***end prologue       Odd_Polynomials

  SUBROUTINE Odd_Polynomials(grid,reg_poly_o,typwt)
  IMPLICIT NONE
  TYPE(regional_grid)                  :: grid
  TYPE(odd_functions),  DIMENSION(:)   :: reg_poly_o
  CHARACTER(LEN=*)                     :: typwt
  REAL*8                               :: one = 1.d0
  REAL*8                               :: dum
  INTEGER                              :: i
!
!
  DO  i=1,nreg
      Write(iout,2) i, npt(i)
!
!        calculate the sector functions and their derivatives for the m odd grid.
!
      ALLOCATE( grid%reg_poly_o(i)%pr(npt(i),npt(i)),                      &
                grid%reg_poly_o(i)%dpr(npt(i),npt(i)),                     &
                grid%reg_poly_o(i)%ddpr(npt(i),npt(i)) )
      Call renorm( grid,i,npt(i) )                               
  END DO
1 FORMAT(/,10x,'Calculate the Regional Basis Functions for Coordinate = ',a4)
2 FORMAT(/,10x,'Region = ',i3,2x,'Number of Polynomials = ', i3)
END SUBROUTINE Odd_Polynomials
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
                                *                                        &
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
           END MODULE Spheroidal_Fedvr_Module
!***********************************************************************
!***********************************************************************
