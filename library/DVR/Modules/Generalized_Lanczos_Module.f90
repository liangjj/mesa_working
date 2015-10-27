!***********************************************************************
! Spheroidal_Fedvr_Module
!**begin prologue     Spheroidal_Fedvrs_Module
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
                           USE dvr_global
                           USE dvr_shared
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                           INTERFACE Ke_Fedvr
             MODULE PROCEDURE Eta_Ke_Even,                             &
                              Eta_Ke_Odd,                              &
                              Xi_Ke_Even,                              &
                              Xi_Ke_Odd
                       END INTERFACE Ke_Fedvr
!
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
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
!***********************************************************************
!***********************************************************************
!deck Eta_KE_Even.f
!***begin prologue     Eta_KE_Even
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            angular KE elements for prolate coordinates.
!*** description       Here we consider the matrix elements
!***                   < f_i | d/dx (1 - x*x ) d/dx | f_j>
!                      = - <  | [ df_i/dx df_j/dx ] ( 1 - x*x ) >
!                      where the last line comes from integrating by parts and dropping the
!                      vanishing surface term.       
!***references

!***routines called
!***end prologue       Eta_KE_Even

  SUBROUTINE Eta_KE_Even(tr,q,wt,f,df,n)
  USE dvr_global ,  ONLY  : iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                              :: n
  REAL*8, DIMENSION(:)                 :: q
  REAL*8, DIMENSION(:)                 :: wt
  REAL*8, DIMENSION(:,:)               :: f
  REAL*8, DIMENSION(:,:)               :: df
  REAL*8, DIMENSION(:,:)               :: tr
  REAL*8, DIMENSION(:)                 :: fac
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: k
  CALL rzero(tr,n*n)
  END DO
  DO i = 1, n
     DO j = 1, i
        DO k = 1, n
           tr(i,j) = tr(i,j) - fac(k) * wt(k) * df(k,i) * df(k,j) 
        END DO
        tr(j,i) = tr(i,j)
     END DO
  END DO
END SUBROUTINE Eta_KE_Even
!***********************************************************************
!***********************************************************************
!deck Eta_Ke_Odd.f
!***begin prologue     Eta_Ke_Odd
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            angular KE elements for prolate coordinates.
!*** description       Here we consider the matrix elements
!***                   < f_i | d/dx (1 - x*x ) d/dx | f_j>
!                      = - <  | [ df_i/dx df_j/dx ] ( 1 - x*x ) >
!                      where the last line comes from integrating by parts and dropping the
!                      vanishing surface term.       
!***references

!***routines called
!***end prologue       Eta_Ke_Odd

  SUBROUTINE Eta_Ke_Odd(tr,q,wt,f,df,n)
  USE dvr_global ,  ONLY  : iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                              :: n
  REAL*8, DIMENSION(:)                 :: q
  REAL*8, DIMENSION(:)                 :: wt
  REAL*8, DIMENSION(:,:)               :: f
  REAL*8, DIMENSION(:,:)               :: df
  REAL*8, DIMENSION(:,:)               :: tr
  REAL*8, DIMENSION(:)                 :: fac
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: k
  CALL rzero(tr,n*n)
  DO  i = 1, n
      DO j = 1, i
         DO k = 1, n
            tr(i,j) = tr(i,j) -  fac(k) * fac(k) * wt(k) * df(k,i) * df(k,j)
         END DO
         tr(i,j) = tr(i,j)  + q(i) * fac(i) * f(i,i) * df(i,j)    &
                            + q(j) * fac(j) * f(j,j) * df(j,i) ) &
         tr(j,i) = tr(i,j)
      END DO
      tr(i,i) = tr(i,i)  - q(i) * q(i) * f(i,i) * f(i,i)
  END DO
END SUBROUTINE Eta_Ke_Odd
!***********************************************************************
!***********************************************************************
!deck Xi_KE_Even.f
!***begin prologue     Xi_KE_Even
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            angular KE elements for prolate coordinates.
!*** description       Here we consider the matrix elements
!***                   < f_i | d/dx (x*x - 1. ) d/dx | f_j>
!                      = - <  | [ df_i/dx df_j/dx ] ( 1 - x*x ) >
!                      where the last line comes from integrating by parts and dropping the
!                      vanishing surface term.       
!***references

!***routines called
!***end prologue       Xi_KE_Even

  SUBROUTINE Xi_KE_Even(tr,q,wt,f,df,n)
  USE dvr_global ,  ONLY  : iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                              :: n
  REAL*8, DIMENSION(:)                 :: q
  REAL*8, DIMENSION(:)                 :: wt
  REAL*8, DIMENSION(:,:)               :: f
  REAL*8, DIMENSION(:,:)               :: df
  REAL*8, DIMENSION(:,:)               :: tr
  REAL*8, DIMENSION(:)                 :: fac
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: k
  CALL rzero(tr,n*n)
  DO  i = 1, n     
      DO j = 1, i
         DO k = 1, n
            tr(i,j) = tr(i,j) - fac(k) * wt(k) * df(k,i) * df(k,j) 
         END DO
         tr(j,i) = tr(i,j)
      END DO
  END DO
END SUBROUTINE Xi_KE_Even
!***********************************************************************
!***********************************************************************
!deck Xi_KE_Odd.f
!***begin prologue     Xi_KE_Odd
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            angular KE elements for prolate coordinates.
!*** description       Here we consider the matrix elements
!***                   < f_i | d/dx (x*x - 1. ) d/dx | f_j>
!                      = - <  | [ df_i/dx df_j/dx ] ( 1 - x*x ) >
!                      where the last line comes from integrating by parts and dropping the
!                      vanishing surface term.       
!***references

!***routines called
!***end prologue       Xi_KE_Odd

  SUBROUTINE Xi_KE_Odd(tr,q,wt,f,df,n)
  USE dvr_global ,  ONLY  : iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                              :: n
  REAL*8, DIMENSION(:)                 :: q
  REAL*8, DIMENSION(:)                 :: wt
  REAL*8, DIMENSION(:,:)               :: f
  REAL*8, DIMENSION(:,:)               :: df
  REAL*8, DIMENSION(:,:)               :: tr
  REAL*8, DIMENSION(:)                 :: fac
  INTEGER                              :: i
  INTEGER                              :: j
  INTEGER                              :: k
  CALL rzero(tr,n*n)
  DO  i = 1, n
      DO j = 1, i
         DO k = 1, n
            tr(i,j) = tr(i,j) -  fac(k) * fac(k) * wt(k) * df(k,i) * df(k,j)
         END DO
         tr(i,j) = tr(i,j)  - q(i) * fac(i) * f(i,i) * df(i,j)    &
                            - q(j) * fac(j) * f(j,j) * df(j,i) )  
         tr(j,i) = tr(i,j)
      END DO
      tr(i,i) = tr(i,i)  - q(i) * q(i) * f(i,i) * f(i,i)
  END DO
END SUBROUTINE Xi_KE_Odd
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
           END MODULE Spheroidal_Fedvr_Module
!***********************************************************************
!***********************************************************************
