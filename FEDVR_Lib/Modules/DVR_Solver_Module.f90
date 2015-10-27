!***********************************************************************
! DVR_Solver_Module
!**begin prologue     DVR_Solver_Module
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            Read in FEDVR Data
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      DVR_Solver_Module
!***********************************************************************
!***********************************************************************
                           MODULE DVR_Solver_Module
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!
                            INTERFACE Solution                       
                       MODULE PROCEDURE Solution_Cartesian,               &
                                        Solution_Spherical,               &
                                        Solution_Spheroidal  
                            END INTERFACE Solution
!
                            INTERFACE Compare                       
                       MODULE PROCEDURE Compare_Cartesian,                &
                                        Compare_Spherical,                &
                                        Compare_Spheroidal  
                            END INTERFACE Compare
!
                            INTERFACE Two_Electron_Q                       
                       MODULE PROCEDURE Two_Electron_Spherical_Q,         &
                                        Two_Electron_Spheroidal_Q  
                            END INTERFACE Two_Electron_Q  
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Diagonalize.f
!***begin prologue     Diagonalize
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            diagonalization of DVR/FEM hamiltonian.
!***
!***references
!***routines called
!***end prologue       Diagonalize
  SUBROUTINE Diagonalize(grid, dvr_mat, val)
  IMPLICIT NONE
  TYPE(coordinates)                      :: grid  
  TYPE (dvr_matrices), DIMENSION(0:)     :: dvr_mat
  INTEGER                                :: info
  INTEGER                                :: i
  INTEGER                                :: j
  INTEGER                                :: val
  INTEGER                                :: count
  CHARACTER(LEN=3)                       :: itoc
!
! Use the packed form of the diagonalization routine
!
  count = 0
  DO i = 1, physical_points
     DO j = 1, i
        count = count + 1
        dvr_mat(0)%lower(count) = dvr_mat(val)%ham(i,j)
     END DO
  END DO
!
  Call dspev('v','u',physical_points,dvr_mat(0)%lower,              &
             dvr_mat(0)%eigenvalues,dvr_mat(0)%eigenvectors,        &
             physical_points,dvr_mat(0)%work,info)
!
!
  title='eigenvalues for angular quantum number = '//itoc(val)
  CALL prntrm(title,dvr_mat(0)%eigenvalues,physical_points,1,physical_points,1,iout)
  IF(prn(10)) THEN
     title='eigenvectors'
     CALL prntrm(title,dvr_mat(0)%eigenvectors,physical_points,physical_points,                       &
                                               physical_points,physical_points,iout)
  END IF
END SUBROUTINE Diagonalize
!***********************************************************************
!***********************************************************************
!deck Solver.f
!***begin prologue     Solver
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Solve the Poisson equation using packed solvers.
!***
!***references
!***routines called
!***end prologue       Solver
  SUBROUTINE Solver(grid, dvr_mat, val)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE (dvr_matrices), DIMENSION(0:)       :: dvr_mat
  REAL(idp)                                :: scale  
  INTEGER                                  :: info
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: val
!
! Using a packed matrix
!
  scale = one / dscale
  count = 0
  DO i = 1, dvr_mat(0)%n_final
     DO j = 1, i
        count = count + 1
        dvr_mat(0)%lower(count) = dvr_mat(val)%ham(i,j)
     END DO
  END DO
  dvr_mat(0)%lower(:) = scale * dvr_mat(0)%lower(:)
!
! Forming the inverse by solving for the unit matrix RHS
! By doing it this way we can compute any density as a RHS
!
  dvr_mat(val)%t(:,:) = zero 
  DO i = 1, dvr_mat(0)%n_final 
     dvr_mat(val)%t(i,i) = one
  END DO
  Call DSPSV('u',dvr_mat(0)%n_final,dvr_mat(0)%n_final,                        &
             dvr_mat(0)%lower,dvr_mat(0)%ipvt,dvr_mat(val)%t,                  &
             dvr_mat(0)%n_final,info)
END SUBROUTINE Solver
!***********************************************************************
!***********************************************************************
!deck Solution_Cartesian
!***begin prologue     Solution_Cartesian
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form the right hand sides of the Poisson equation
!***                   for some special test cases basically in cartesian coordinates.
!***                   
!***                   
!***                   
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Solution_Cartesian
!
  SUBROUTINE Solution_Cartesian(grid,dvr_mat,val,name_cartesian)
  IMPLICIT NONE
  TYPE(coordinates)                   :: grid  
  TYPE(cartesian)                     :: name_cartesian
  TYPE (dvr_matrices), DIMENSION(0:)  :: dvr_mat
  INTEGER                             :: i
  INTEGER                             :: val
  IF (dvr_mat(0)%type_inhomo == 'one') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = one 
      END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'linear') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = grid%grid_points(1:dvr_mat(0)%n_final)  
      END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'exponential') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = exp(-grid%grid_points(1:dvr_mat(0)%n_final))
      END DO
  ELSE IF (dvr_mat(0)%type_inhomo == 'harmonic') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = exp( -grid%grid_points(1:dvr_mat(0)%n_final)     &
                                                                           *                                       &
                                                                  grid%grid_points(1:dvr_mat(0)%n_final)) 
      END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'coulomb') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = one / grid%grid_points(1:dvr_mat(0)%n_final) 
      END DO
  END IF
  DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
     dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) &
                                                                      *                                              &
                                                              sqrt( grid%grid_weights(1:dvr_mat(0)%n_final) )
  END DO
!
! Use exact solution as temporary storage.
!
  dvr_mat(val)%exact_solution = zero
  call ebcx ( dvr_mat(val)%exact_solution, physical_points,                                                        &
              dvr_mat(val)%t, physical_points,                                                                     &
              dvr_mat(val)%computed_solution, physical_points,                                                     &
              dvr_mat(0)%n_final, dvr_mat(0)%n_final, dvr_mat(0)%number_of_right_hand_sides )
  dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,:) = dvr_mat(val)%exact_solution(1:dvr_mat(0)%n_final,:) 
END SUBROUTINE Solution_Cartesian
!***********************************************************************
!***********************************************************************
!deck Solution_Spherical
!***begin prologue     Solution_Spherical
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form the right hand sides of the Poisson equation
!***                   for the spherical coordinate system.  In this case
!***                   we need the full inverse.
!***                   
!***                   
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Solution_Spherical
!
  SUBROUTINE Solution_Spherical(grid,dvr_mat,val,name_spherical)
  IMPLICIT NONE
  TYPE(coordinates)                   :: grid  
  TYPE(spherical)                     :: name_spherical
  TYPE (dvr_matrices), DIMENSION(0:)  :: dvr_mat
  INTEGER                             :: i
  INTEGER                             :: val
  INTEGER                             :: l_factor
  REAL(idp), dimension(3,3)           :: a, b, c
!
  l_factor = - ( val + val + int_one )
  IF (dvr_mat(0)%type_inhomo == 'one') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = one           
      END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'linear') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = grid%grid_points(1:dvr_mat(0)%n_final)           
      END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'exponential') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = exp(-grid%grid_points(1:dvr_mat(0)%n_final))     
      END DO
  ELSE IF (dvr_mat(0)%type_inhomo == 'harmonic') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = exp( -grid%grid_points(1:dvr_mat(0)%n_final)     &
                                                                           *                                       &
                                                                  grid%grid_points(1:dvr_mat(0)%n_final)) 
      END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'coulomb') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = one / grid%grid_points(1:dvr_mat(0)%n_final)     
      END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'basis') THEN
      dvr_mat(val)%computed_solution(:,:) = zero 
      DO i = 1, dvr_mat(0)%n_final 
         dvr_mat(val)%computed_solution(i,i) =     one        
      END DO
  END IF
  IF(dvr_mat(0)%type_inhomo /= 'basis') THEN
     DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
        dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) &
                                                                             *                                          &
                                                                   sqrt( grid%grid_weights(1:dvr_mat(0)%n_final) )      &
                                                                             *                                          &
                                                                   grid%grid_points(1:dvr_mat(0)%n_final)       
     END DO
  ELSE
     DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
        dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) &
                                                                             *                                          &
                                                                         grid%grid_points(1:dvr_mat(0)%n_final)         &
                                                                             /                                          &
                                                                         sqrt( grid%grid_weights(1:dvr_mat(0)%n_final) ) 
     END DO
     dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,:) = l_factor *                                                &
                                                              dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,:) 
  END IF
!
! Use exact solution as temporary storage.
!
  call ebc  ( dvr_mat(val)%exact_solution,                                         &
              dvr_mat(val)%t,                                                      &
              dvr_mat(val)%computed_solution,                                      &
              dvr_mat(0)%n_final, dvr_mat(0)%n_final, dvr_mat(0)%number_of_right_hand_sides )
  dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,:) = dvr_mat(val)%exact_solution(1:dvr_mat(0)%n_final,:)
END SUBROUTINE Solution_Spherical
!***********************************************************************
!***********************************************************************
!deck Solution_Spheroidal
!***begin prologue     Solution_Spheroidal
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form the right hand sides of the Poisson equation
!***                   for some special test cases.
!***                   
!***                   
!***                   
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Solution_Spheroidal
!
  SUBROUTINE Solution_Spheroidal(grid,dvr_mat,val,name_spheroidal)
  IMPLICIT NONE
  TYPE(coordinates)                   :: grid  
  TYPE(spheroidal)                    :: name_spheroidal  
  TYPE (dvr_matrices), DIMENSION(0:)  :: dvr_mat
  INTEGER                             :: i
  INTEGER                             :: val
  IF (dvr_mat(0)%type_inhomo == 'one') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = one 
      END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'linear') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = grid%grid_points(1:dvr_mat(0)%n_final)  
      END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'exponential') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = exp(-grid%grid_points(1:dvr_mat(0)%n_final))
      END DO
  ELSE IF (dvr_mat(0)%type_inhomo == 'harmonic') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = exp( -grid%grid_points(1:dvr_mat(0)%n_final)   &
                                                                            *                                    &
                                                                  grid%grid_points(1:dvr_mat(0)%n_final)) 
      END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'coulomb') THEN
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
         dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = one / grid%grid_points(1:dvr_mat(0)%n_final) 
      END DO
  END IF
  DO i = 1, dvr_mat(0)%number_of_right_hand_sides 
     dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) = dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i) &
                                                                                *                                    &
                                                              sqrt( grid%grid_weights(1:dvr_mat(0)%n_final) )
  END DO
END SUBROUTINE Solution_Spheroidal

!***********************************************************************
!***********************************************************************
!deck Compare_Cartesian.f
!***begin prologue     Compare_Cartesian
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Compare computed and exact solutions to the Poisson
!***                   equation for some specific cases.
!***
!***references
!***routines called
!***end prologue       Compare_Cartesian
  SUBROUTINE Compare_Cartesian(grid, dvr_mat,val, name_cartesian)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE (dvr_matrices), DIMENSION(0:)       :: dvr_mat
  TYPE(cartesian)                          :: name_cartesian
  REAL(idp)                                :: added  
  INTEGER                                  :: info
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: val
  CHARACTER(LEN=3)                         :: itoc
  CHARACTER(LEN=80)                        :: title
  DO i = 1, dvr_mat(0)%number_of_right_hand_sides
     dvr_mat(val)%t(1:dvr_mat(0)%n_final,i) = dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i)   &
                                                             /                                         &
                                              sqrt( grid%grid_weights(1:dvr_mat(0)%n_final) )  
  END DO
!
!           Some tests for cases on (0,1) where the solution vanishes at the left end
!           and the solution either vanishes or is one at the right end.
!
  IF (dvr_mat(0)%type_inhomo == 'one') THEN
!
!         Particular solution which vanishes at right end x *  ( x - 1 ) / 2
!
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides
         dvr_mat(val)%exact_solution(:,i) =  grid%grid_points(:) *                  &
                                             ( grid%grid_points(:) - one ) * half
      END DO
!
  ELSE IF(dvr_mat(0)%type_inhomo == 'linear') THEN
!
!         Solution is x *  ( x*x - 1 ) / 6
!
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides
         dvr_mat(val)%exact_solution(:,i) =  grid%grid_points(:) *                  &
                                           ( grid%grid_points(:) *                  &
                                             grid%grid_points(:) - one ) * sixth
      END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'exponential') THEN
!
!         Solution is exp(-x) + x *  ( 1 -exp(-1) ) - 1  
!
      added = exp(-one)
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides
         dvr_mat(val)%exact_solution(:,i) =  exp(-grid%grid_points(:))              &
                                                     +                              &
                                             grid%grid_points(:)                    &
                                                     *                              &
                                             ( one - added ) - one                    
      END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'coulomb') THEN      
!
!         Solution is x * ( ln (x) - ln(1) )
!
      added = log(one)
      DO i = 1, dvr_mat(0)%number_of_right_hand_sides
         dvr_mat(val)%exact_solution(:,i) =  grid%grid_points(:)                    &
                                                  *                                 &
                                            ( log ( grid%grid_points(:) ) - added )
      END DO
  ELSE
      Call lnkerr('not valid')
  END IF
!
  IF (.not.drop(2)) THEN
      DO i= 1, dvr_mat(0)%n_final
         dvr_mat(val)%t(i,:) = dvr_mat(val)%t(i,:) + grid%grid_points(i)       
         dvr_mat(val)%exact_solution(i,:)    = dvr_mat(val)%exact_solution(i,:)    + grid%grid_points(i)       
      END DO
      dvr_mat(val)%t(physical_points,:) = R_max             
      dvr_mat(val)%exact_solution(physical_points,:)    = R_max             
  END IF
  title='Numerical Solution'
  Call Prntfm(title,dvr_mat(val)%t,physical_points,                                                &
                    dvr_mat(0)%number_of_right_hand_sides, physical_points,                        &
                    dvr_mat(0)%number_of_right_hand_sides,iout,'e')
  title='Exact Solution'
  Call Prntfm(title,dvr_mat(val)%exact_solution,physical_points,                                   &
                    dvr_mat(0)%number_of_right_hand_sides, physical_points,                        &
                    dvr_mat(0)%number_of_right_hand_sides,iout,'e')
END SUBROUTINE Compare_Cartesian
!***********************************************************************
!***********************************************************************
!deck Compare_Spherical
!***begin prologue     Compare_Spherical
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Compare computed and exact solutions to the Poisson
!***                   equation for some specific cases.
!***
!***references
!***routines called
!***end prologue       Compare_Spherical
  SUBROUTINE Compare_Spherical(grid, dvr_mat, val, name_spherical)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE (dvr_matrices), DIMENSION(0:)       :: dvr_mat
  TYPE(spherical)                          :: name_spherical
  REAL(idp)                                :: scale  
  INTEGER                                  :: info
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: val
  CHARACTER(LEN=3)                         :: itoc
  CHARACTER(LEN=80)                        :: title
!
! In the spherical case the reason for dividing by the grid points is because we have actually solved the
! equations where the volume element r**2 scaled the original overall set of equations.  By transforming this
! out of the problem, the inverse solved for actually has an additional factor of r in its definition.  This
! transforms it back by dividing out that factor.
!
  DO i = 1, dvr_mat(0)%number_of_right_hand_sides
     dvr_mat(val)%t(1:dvr_mat(0)%n_final,i) = dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i)   &
                                                             /                                         &
                                                    ( sqrt( grid%grid_weights(1:dvr_mat(0)%n_final) )  &
                                                                      *                                &
                                                            grid%grid_points(1:dvr_mat(0)%n_final) ) 
  END DO
!
!         Some tests for cases where solution is one at r = R
!
!
  IF (dvr_mat(0)%type_inhomo == 'one') THEN
!
!            Solution is (r * r - R * R) / 6 
!
     DO i = 1, dvr_mat(0)%number_of_right_hand_sides
        dvr_mat(val)%exact_solution(:,i) =  ( grid%grid_points(:)                  &
                                                      *                            &
                                              grid%grid_points(:)                  &
                                                      -                            &
                                                     R_max                         &
                                                      *                            &
                                                     R_max )                       &
                                                             * sixth 
     END DO
!
  ELSE IF(dvr_mat(0)%type_inhomo == 'linear') THEN
!
!            Solution is (r * r * r - R * R * R) / 12 
!
     DO i = 1, dvr_mat(0)%number_of_right_hand_sides
        dvr_mat(val)%exact_solution(:,i) =  ( grid%grid_points(:)                  &
                                                       *                           &
                                              grid%grid_points(:)                  &
                                                       *                           &
                                              grid%grid_points(:)                  &
                                                       -                           &
                                                       R_max                       &
                                                        *                          &
                                                       R_max                       &
                                                        *                          &
                                                       R_max ) * sixth * half
     END DO
  ELSE IF(dvr_mat(0)%type_inhomo == 'coulomb') THEN      
!
!            Solution is (r - R) / 2 
!
     DO i = 1, dvr_mat(0)%number_of_right_hand_sides
        dvr_mat(val)%exact_solution(:,i) =  ( grid%grid_points(:) - R_max )       &
                                                  *                               &
                                                 half
     END DO
  ELSE
!
     Call lnkerr('not valid')
!
  END IF
!
  IF (.not.drop(2)) THEN
      DO i= 1, dvr_mat(0)%n_final
         dvr_mat(val)%t(i,:) = dvr_mat(val)%t(i,:) + grid%grid_points(i)       
         dvr_mat(val)%exact_solution(i,:)    = dvr_mat(val)%exact_solution(i,:)    + grid%grid_points(i)       
      END DO
      dvr_mat(val)%t(physical_points,:) = R_max             
      dvr_mat(val)%exact_solution(physical_points,:)    = R_max             
  END IF
  title='Numerical Solution for angular quantum number = '//itoc(val)
  Call Prntfm(title,dvr_mat(val)%t,physical_points,                                                &
                    dvr_mat(0)%number_of_right_hand_sides, physical_points,                      &
                    dvr_mat(0)%number_of_right_hand_sides,iout,'e')
  title='Exact Solution for angular quantum number = '//itoc(val)
  Call Prntfm(title,dvr_mat(val)%exact_solution,physical_points,                                   &
                    dvr_mat(0)%number_of_right_hand_sides, physical_points,                      &
                    dvr_mat(0)%number_of_right_hand_sides,iout,'e')
END SUBROUTINE Compare_Spherical
!***********************************************************************
!***********************************************************************
!deck Compare_Spheroidal
!***begin prologue     Compare_Spheroidal
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Compare computed and exact solutions to the Poisson
!***                   equation for some specific cases.
!***
!***references
!***routines called
!***end prologue       Compare_Spheroidal
  SUBROUTINE Compare_Spheroidal(grid, dvr_mat,val, name_spheroidal)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE (dvr_matrices), DIMENSION(0:)       :: dvr_mat
  TYPE(spheroidal)                         :: name_spheroidal
  REAL(idp)                                :: scale  
  INTEGER                                  :: info
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: val
  CHARACTER(LEN=3)                         :: itoc
  CHARACTER(LEN=80)                        :: title
  DO i = 1, dvr_mat(0)%number_of_right_hand_sides
     dvr_mat(val)%rhs(1:dvr_mat(0)%n_final,i) = dvr_mat(val)%computed_solution(1:dvr_mat(0)%n_final,i)   &
                                                         /                                           &
                                              sqrt( grid%grid_weights(1:dvr_mat(0)%n_final) )
  END DO
END SUBROUTINE Compare_Spheroidal
!***********************************************************************
!***********************************************************************
!deck Two_Electron_Spherical_Q
!***begin prologue     Two_Electron_Spherical_Q
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Compare computed and exact solutions to the Poisson
!***                   equation for some specific cases.
!***
!***references
!***routines called
!***end prologue       Two_Electron_Spherical_Q
  SUBROUTINE Two_Electron_Spherical_Q (grid, dvr_mat, val, name_spherical)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE (dvr_matrices), DIMENSION(0:)       :: dvr_mat
  TYPE(spherical)                          :: name_spherical
  REAL(idp)                                :: R_factor  
  REAL(idp)                                :: fac_1
  REAL(idp)                                :: fac_2  
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: val
  INTEGER                                  :: l_factor  
  CHARACTER(LEN=3)                         :: itoc
  CHARACTER(LEN=80)                        :: title
!
! Since we are interested in the solution where the right hand side is the unit matrix, we need to make
! sure the r**2 in the volume is handled properly.  This ensure that to be the case.
  DO i = 1, physical_points
     DO j = 1, physical_points
        dvr_mat(val)%t(i,j) = dvr_mat(val)%t(i,j) / ( grid%grid_points(i) *  grid%grid_points(j) )  
     END DO
  END DO
!
!       Q^L(i,j) = ( r_i )^(L+2)*( r_j )^L /( r_n)^(2L+1) - (2L+1) * ( r_i )^2 * ( T_ji )^-1 /sqrt(w_i*w_j}
!
  l_factor = ( val + val + int_one )
  R_factor = one / R_max**l_factor
  DO i = 1, physical_points
     fac_1 = ( grid%grid_points(i) ) ** ( val + int_two ) * R_factor 
     fac_2 = l_factor * ( grid%grid_points(i) ) ** int_two
     DO j = 1, physical_points
        dvr_mat(val)%Q(i,j) = ( grid%grid_points(j) ** val ) * fac_1 - fac_2 * dvr_mat(val)%t(j,i)     &
                                                       /                                               &
                                ( sqrt ( grid%grid_weights(i) * grid%grid_weights(j) ) )
     END DO
  END DO
  title='Q Function for Two Electron Integrals_'//itoc(val)
  Call Prntfm(title,dvr_mat(val)%Q,physical_points, physical_points, physical_points,                  &
              physical_points, iout,'e')
END SUBROUTINE Two_Electron_Spherical_Q
!***********************************************************************
!***********************************************************************
!deck Two_Electron_Spheroidal_Q
!***begin prologue     Two_Electron_Spheroidal_Q
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Compare computed and exact solutions to the Poisson
!***                   equation for some specific cases.
!***
!***references
!***routines called
!***end prologue       Two_Electron_Spheroidal_Q
  SUBROUTINE Two_Electron_Spheroidal_Q(grid, dvr_mat, val, name_spheroidal)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE (dvr_matrices), DIMENSION(0:)       :: dvr_mat
  TYPE(spheroidal)                         :: name_spheroidal
  REAL(idp)                                :: scale  
  INTEGER                                  :: info
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: val
  CHARACTER(LEN=3)                         :: itoc
  CHARACTER(LEN=80)                        :: title
!
END SUBROUTINE Two_Electron_Spheroidal_Q
!***********************************************************************
!***********************************************************************
           END MODULE DVR_Solver_Module
!***********************************************************************
!***********************************************************************
