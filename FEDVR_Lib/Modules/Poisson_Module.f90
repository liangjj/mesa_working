!***********************************************************************
! Poisson_Module
!**begin prologue     Poisson_Module
!**date written       090219   (yymmdd)
!**revision date               (yymmdd)
!**keywords           DVR, FEDVR
!**
!**author             schneider, b. i.(nsf)
!**source             DVR Library
!**purpose            Driver module to calculate the FEDVR functions and matrix elements
!***                  in a FEDVR basis
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Poisson_Module
!***********************************************************************
!***********************************************************************
                           MODULE Poisson_Module
                           USE Data_Module
                           USE FEDVR_Shared
                           USE FEDVR_Derived_Types
                           USE Legendre
                     INTEGER                      :: len_1
!***********************************************************************
!                           Explicit Interfaces
!***********************************************************************
!
                            INTERFACE Right_Hand_Side                       
                       MODULE PROCEDURE Right_Hand_Side_Cartesian,        &
                                        Right_Hand_Side_Spherical,        &
                                        Right_Hand_Side_Spheroidal  
                            END INTERFACE Right_Hand_Side
!
                            INTERFACE Final_Solution                       
                       MODULE PROCEDURE Final_Solution_Cartesian,         &
                                        Final_Solution_Spherical,         &
                                        Final_Solution_Spheroidal  
                            END INTERFACE Final_Solution
!
                            INTERFACE Exact                       
                       MODULE PROCEDURE Exact_Cartesian,                  &
                                        Exact_Spherical,                  &
                                        Exact_Spheroidal  
                            END INTERFACE Exact
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Poisson_Equation
!***begin prologue     Poisson_Equation
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Using the Hamiltonian solve the Poisson equation.
!***                   
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Poisson_Equation
!
  SUBROUTINE Poisson_Equation(grid)
  IMPLICIT NONE
  TYPE (coordinates)                             :: grid
  INTEGER                                        :: i
  INTEGER                                        :: lm
  INTEGER                                        :: l_factor
  REAL(idp)                                      :: fpkey
  LOGICAL                                        :: dollar
  LOGICAL                                        :: logkey
  INTEGER                                        :: intkey
  CHARACTER(LEN=80)                              :: chrkey
  CHARACTER(LEN=3)                               :: itoc
!
!
!
  IF ( dollar('$poisson_'//grid%label(1:len),card,cpass,inp) ) THEN
       number_of_right_hand_sides=intkey(card,'number_of_right_hand_sides',            &
                                               physical_points,' ')
       type_inhomo=chrkey(card,'inhomogeneity','linear',' ') 
       last_value = fpkey(card,'last_value',one,' ')
  END IF
!
!      The number of physical points may or may not contain the last boundary point.
!      If the eigenvalue problem solved was for a solution which is zero at the boundary
!      point, then it was dropped and the Hamiltonian matrix does not contain that
!      point.  Under this circumstance, the number of global points is larger by one
!      than the number of physical points.

!
  n_final=physical_points
  n_total = physical_points
  IF(.not.drop(2)) THEN
      n_final = physical_points - int_one
  ELSE
      n_total = n_total + int_one
  END IF
  ALLOCATE( dvr_mat(0)%ipvt( n_total ),                                          &
            dvr_mat(0)%lower( n_total*(n_total+1)/2 ),                           &
            dvr_mat(0)%RHS( n_total,                                             &
                            number_of_right_hand_sides ),                        &        
            dvr_mat(0)%Exact_Solution( n_total,                                  &
                                       number_of_right_hand_sides) )  
!  
  DO lm = 0, size 
!
!    We compute the RHS and then the exact solution for a few selected cases      
!    just as test.  Then we solve the linear equations and compare the two.
!    The numerical solutions are computed for the specific case where the solution
!    vanishes at the right hand boundary.  The a solution of the homogeneous equation
!    is added to satisfy the bounday contitions.
!
     dvr_mat(0)%RHS(:,:) = zero
     IF (keyword == 'cartesian') THEN
         Call Right_Hand_Side(grid,grid%name_cartesian)
         Call Solver(grid,lm)
         Call Final_Solution(grid,grid%name_cartesian,lm)
         title='Numerical Solution to Linear Equations'
         Call Prntfm(title,dvr_mat(0)%RHS,n_total,                             &
                     number_of_right_hand_sides, n_total,                      &
                     number_of_right_hand_sides,iout,'e')
         Call Exact(grid,grid%name_cartesian)
         title='Exact Solution to Linear Equations'
         Call Prntfm(title,dvr_mat(0)%Exact_Solution,n_total,                  &
                     number_of_right_hand_sides, n_total,                      &
                     number_of_right_hand_sides,iout,'e')
     ELSE IF(keyword == 'spherical') THEN
         l_factor = - ( lm + lm + int_one )
         Call Right_Hand_Side(grid,grid%name_spherical)
         Call Solver(grid,lm)
         Call Final_Solution(grid,grid%name_spherical,lm)
         title='Numerical Solution to Linear Equations for '//                 &
               'Angular Quantum Number = '//itoc(lm)
         Call Prntfm(title,dvr_mat(0)%RHS,n_total,                             &
                     number_of_right_hand_sides, n_total,                      &
                     number_of_right_hand_sides,iout,'e')
         Call Exact(grid,grid%name_spherical)
         title='Exact Solution to Linear Equations for '//                     &
               'Angular Quantum Number = '//itoc(lm)
         Call Prntfm(title,dvr_mat(0)%Exact_Solution,n_total,                  &
                     number_of_right_hand_sides, n_total,                      &
                     number_of_right_hand_sides,iout,'e')
     ELSE IF(keyword == 'spheroidal') THEN
         Call Right_Hand_Side(grid,grid%name_spheroidal)
         Call Solver(grid,lm)
         Call Final_Solution(grid,grid%name_spheroidal,lm)
         title='Numerical Solution to Linear Equations for '//                 &
               'Angular Quantum Number = '//itoc(lm)
         Call Prntfm(title,dvr_mat(0)%RHS,n_total,                             &
                     number_of_right_hand_sides, n_total,                      &
                     number_of_right_hand_sides,iout,'e')
         Call Exact(grid,grid%name_spheroidal)
         title='Exact Solution to Linear Equations for '//                     &
               'Angular Quantum Number = '//itoc(lm)
         Call Prntfm(title,dvr_mat(0)%Exact_Solution,n_total,                  &
                     number_of_right_hand_sides, n_total,                      &
                     number_of_right_hand_sides,iout,'e')
     END IF
     IF (.not.two_electron) THEN
         DEALLOCATE(dvr_mat(lm)%tr )
     END IF
  END DO
  Call IOsys('rewind all on '//FEDVR_File//' read-and-write',0,0,0,' ')
  DEALLOCATE( dvr_mat(0)%ipvt,                                                 &
              dvr_mat(0)%lower,                                                &
              dvr_mat(0)%RHS,                                                  &           
              dvr_mat(0)%Exact_Solution )  
!
!
END SUBROUTINE Poisson_Equation
!***********************************************************************
!***********************************************************************
!deck Right_Hand_Side_Cartesian
!***begin prologue     Right_Hand_Side_Cartesian
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
!***end prologue       Right_Hand_Side_Cartesian
!
  SUBROUTINE Right_Hand_Side_Cartesian(grid,name_cartesian)
  IMPLICIT NONE
  TYPE(coordinates)                   :: grid  
  TYPE(cartesian)                     :: name_cartesian
  INTEGER                             :: i

!
!         The weight from the integration is multiplied in
!
  IF (type_inhomo == 'one') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = sqrt( grid%grid_weights(1:n_final) )
      END DO
  ELSE IF(type_inhomo == 'linear') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = grid%grid_points(1:n_final)                            &
                                                         *                                    &
                                               sqrt( grid%grid_weights(1:n_final) )
      END DO
  ELSE IF(type_inhomo == 'exponential') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = exp(-grid%grid_points(1:n_final))                      &
                                                            *                                 &
                                               sqrt( grid%grid_weights(1:n_final) )
      END DO
  ELSE IF (type_inhomo == 'harmonic') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = exp( -grid%grid_points(1:n_final)                      &
                                                            *                                 &
                                               grid%grid_points(1:n_final))                   &
                                                            *                                 &
                                               sqrt( grid%grid_weights(1:n_final) )
     END DO
  ELSE IF(type_inhomo == 'coulomb') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = sqrt( grid%grid_weights(1:n_final) )                   &
                                                  /                                           &
                                               grid%grid_points(1:n_final) 
      END DO
  END IF
END SUBROUTINE Right_Hand_Side_Cartesian
!***********************************************************************
!***********************************************************************
!deck Right_Hand_Side_Spherical
!***begin prologue     Right_Hand_Side_Spherical
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
!***end prologue       Right_Hand_Side_Spherical
!
  SUBROUTINE Right_Hand_Side_Spherical(grid,name_spherical)
  IMPLICIT NONE
  TYPE(coordinates)                   :: grid  
  TYPE(spherical)                     :: name_spherical
  INTEGER                             :: i
  REAL(idp), dimension(3,3)           :: a, b, c
!
!          It is important to remember that the Laplacian has been scaled, so that 
!          the right hand side is actually r time the right hand side. The weights
!          have been multiplied in.
!

  IF (type_inhomo == 'one') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = sqrt( grid%grid_weights(1:n_final) )                    &
                                                            *                                  &
                                               grid%grid_points(1:n_final)       
      END DO
  ELSE IF(type_inhomo == 'linear') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = sqrt( grid%grid_weights(1:n_final) )                    &
                                                           *                                   &
                                               grid%grid_points(1:n_final)                     &
                                                           *                                   &
                                               grid%grid_points(1:n_final)           
      END DO
  ELSE IF(type_inhomo == 'exponential') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = sqrt( grid%grid_weights(1:n_final) )                    &
                                                           *                                   &
                                               grid%grid_points(1:n_final)                     &
                                                           *                                   &
                                               exp(-grid%grid_points(1:n_final))     
      END DO
  ELSE IF (type_inhomo == 'harmonic') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = sqrt( grid%grid_weights(1:n_final) )                    &
                                                           *                                   &
                                               grid%grid_points(1:n_final)                     &
                                                           *                                   &
                                               exp( -grid%grid_points(1:n_final)               &
                                                                  *                            &
                                                     grid%grid_points(1:n_final)) 
      END DO
  ELSE IF(type_inhomo == 'coulomb') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = sqrt( grid%grid_weights(1:n_final) )     
      END DO
  ELSE IF(type_inhomo == 'basis') THEN
!
!           This is the case where there are two DVR functions 
!           as the right hand side
!
     dvr_mat(0)%RHS(:,:) = zero 
     DO i = 1, number_of_right_hand_sides 
        dvr_mat(0)%RHS(1:n_final,i) = l_factor                                                 &
                                                *                                              &
                                              grid%grid_points(1:n_final)                      &
                                               /                                               &
                                                sqrt( grid%grid_weights(1:n_final) ) 
     END DO
  END IF
END SUBROUTINE Right_Hand_Side_Spherical
!***********************************************************************
!***********************************************************************
!deck Right_Hand_Side_Spheroidal
!***begin prologue     Right_Hand_Side_Spheroidal
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
!***end prologue       Right_Hand_Side_Spheroidal
!
  SUBROUTINE Right_Hand_Side_Spheroidal(grid,name_spheroidal)
  IMPLICIT NONE
  TYPE(coordinates)                   :: grid  
  TYPE(spheroidal)                    :: name_spheroidal  
  INTEGER                             :: i
  IF (type_inhomo == 'one') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = sqrt( grid%grid_weights(1:n_final) ) 
      END DO
  ELSE IF(type_inhomo == 'linear') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = sqrt( grid%grid_weights(1:n_final) )                  &
                                                       *                                     &
                                               grid%grid_points(1:n_final)  
      END DO
  ELSE IF(type_inhomo == 'exponential') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = sqrt( grid%grid_weights(1:n_final) )                  &
                                                       *                                     &
                                               exp(-grid%grid_points(1:n_final))
      END DO
  ELSE IF (type_inhomo == 'harmonic') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = sqrt( grid%grid_weights(1:n_final) )                  &
                                                        *                                    &
                                               exp( -grid%grid_points(1:n_final)             &
                                                        *                                    &
                                               grid%grid_points(1:n_final)) 
      END DO
  ELSE IF(type_inhomo == 'coulomb') THEN
      DO i = 1, number_of_right_hand_sides 
         dvr_mat(0)%RHS(1:n_final,i) = sqrt( grid%grid_weights(1:n_final) )                  &
                                                       /                                     &
                                               grid%grid_points(1:n_final) 
      END DO
  END IF
END SUBROUTINE Right_Hand_Side_Spheroidal
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
  SUBROUTINE Solver(grid,lm)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  REAL(idp)                                :: scale  
  INTEGER                                  :: info
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: lm
  CHARACTER(LEN=3)                         :: itoc
!
! Using a packed matrix
!
  count = 0
  DO i = 1, n_final
     DO j = 1, i
        count = count + 1
        dvr_mat(0)%lower(count) = dvr_mat(lm)%tr(i,j)
     END DO
  END DO
!
! Solve the linear equations and then write out the factorization
! for later use if required.  Remember, the solution does not give
! the solution at a point.  It gives the spectral coefficient.
!
  Call DSPTRF ('u',n_final,dvr_mat(0)%lower,dvr_mat(0)%ipvt,info)
!
  Call IOsys('write real "Factored Matrix l = '//itoc(lm)//' " to '//FEDVR_File, &
              count,dvr_mat(0)%lower,0,' ')
  Call IOsys('write integer "Pivot Array l = '//itoc(lm)//' " to '//FEDVR_File, &
              n_final,dvr_mat(0)%ipvt,0,' ')
!
  Call DSPTRS('u',n_final,number_of_right_hand_sides,dvr_mat(0)%lower,         &
               dvr_mat(0)%ipvt,dvr_mat(0)%RHS,n_total,info)

END SUBROUTINE Solver
!***********************************************************************
!***********************************************************************
!deck Final_Solution_Cartesian.f
!***begin prologue     Final_Solution_Cartesian
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Solve the Poisson equation using packed solvers.
!***
!***references
!***routines called
!***end prologue       Final_Solution_Cartesian
  SUBROUTINE Final_Solution_Cartesian(grid, name_cartesian, lm)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE(cartesian)                          :: name_cartesian
  INTEGER                                  :: i
  INTEGER                                  :: lm
  CHARACTER(LEN=3)                         :: itoc
!
!     First convert the expansion coefficients from the spectral to the
!     coordinate basis. Then add the homogeneous solution for the BC
!
  DO i = 1, number_of_right_hand_sides
     dvr_mat(0)%RHS(1:n_final,i) = dvr_mat(0)%RHS(1:n_final,i)          &
                                            /                           &
                                   sqrt( grid%grid_weights(1:n_final) )
  END DO
  DO i= 1, n_final
     dvr_mat(0)%RHS(i,:)  = dvr_mat(0)%RHS(i,:) + grid%grid_points(i)       
  END DO
  dvr_mat(0)%RHS(n_total,:) = one             
  Call IOsys('write real "Cartesian Poisson Solution Matrix L = '       &
             //itoc(lm)//'" to '//FEDVR_File,                             &
             n_total*number_of_right_hand_sides,dvr_mat(0)%RHS,0,' ')
END SUBROUTINE Final_Solution_Cartesian
!***********************************************************************
!***********************************************************************
!deck Final_Solution_Spherical.f
!***begin prologue     Final_Solution_Spherical
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Solve the Poisson equation using packed solvers.
!***
!***references
!***routines called
!***end prologue       Final_Solution_Spherical
  SUBROUTINE Final_Solution_Spherical(grid, name_spherical, lm)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE(spherical)                          :: name_spherical
  INTEGER                                  :: i
  INTEGER                                  :: lm
  CHARACTER(LEN=3)                         :: itoc
!
!     First convert the expansion coefficients from the spectral to the
!     coordinate basis. Then add the homogeneous solution for the BC
!     The division by the grid point is because we actually used the scaled
!     Laplacian to remove the r**2 factor from the volume element.
!
  DO i =1, number_of_right_hand_sides
     dvr_mat(0)%RHS(1:n_final,i) = dvr_mat(0)%RHS(1:n_final,i) / ( sqrt( grid%grid_weights(:) )     &
                                                         *                                          &
                                                   grid%grid_points(:) )
  END DO
  IF(type_inhomo == 'basis') THEN
     DO i = 1, number_of_right_hand_sides
        last_value = grid%grid_points(i) ** ( lm + int_two) / ( R_max ** ( lm + lm + int_one ) )            
        dvr_mat(0)%RHS(:,i) = last_value * grid%grid_points(:) ** lm + dvr_mat(0)%RHS(:,i)           

     END DO
  ELSE
     dvr_mat(0)%RHS(:,:) = dvr_mat(0)%RHS(:,:) + last_value
  END IF
  Call IOsys('write real "Spherical Poisson Solution Matrix L = '//itoc(lm)//'" to '//FEDVR_File,     &
              n_total*number_of_right_hand_sides,dvr_mat(0)%RHS,0,' ')
END SUBROUTINE Final_Solution_Spherical
!***********************************************************************
!***********************************************************************
!deck Final_Solution_Spheroidal.f
!***begin prologue     Final_Solution_Spheroidal
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Solve the Poisson equation using packed solvers.
!***
!***references
!***routines called
!***end prologue       Final_Solution_Spheroidal
  SUBROUTINE Final_Solution_Spheroidal(grid, name_spheroidal, lm)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE(spheroidal)                         :: name_spheroidal
  INTEGER                                  :: i
  INTEGER                                  :: lm
  CHARACTER(LEN=3)                         :: itoc
!
! Add the homogeneous solution for the BC
!
  DO i =1, number_of_right_hand_sides
     dvr_mat(0)%RHS(:,i) = dvr_mat(0)%RHS(:,i) / sqrt( grid%grid_weights(:) )     
  END DO
  dvr_mat(0)%RHS(n_total,:) = last_value
  Call IOsys('write real "Spheroidal Poisson Solution Matrix L = '      &
             //itoc(lm)//'" to '//FEDVR_File,                             &
             n_total*number_of_right_hand_sides,dvr_mat(0)%RHS,0,' ')
END SUBROUTINE Final_Solution_Spheroidal
!***********************************************************************
!***********************************************************************
!deck Exact_Cartesian.f
!***begin prologue     Exact_Cartesian
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Exact solution to the Poisson
!***                   equation for some specific cases.
!***
!***references
!***routines called
!***end prologue       Exact_Cartesian
  SUBROUTINE Exact_Cartesian(grid, name_cartesian)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE(cartesian)                          :: name_cartesian
  REAL(idp)                                :: added  
  INTEGER                                  :: info
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  CHARACTER(LEN=3)                         :: itoc
  CHARACTER(LEN=80)                        :: title
!
!           Some tests for cases on (0,1) where the solution vanishes at the left end
!           and the solution either vanishes or is one at the right end.
!
  IF (type_inhomo == 'one') THEN
!
!         Particular solution which vanishes at right end x *  ( x - 1 ) / 2
!
      DO i = 1, number_of_right_hand_sides
         dvr_mat(0)%Exact_Solution(1:n_final,i) =  grid%grid_points(:) *                            &
                                                 ( grid%grid_points(:) - one ) * half
      END DO
!
  ELSE IF(type_inhomo == 'linear') THEN
!
!         Solution is x *  ( x*x - 1 ) / 6
!
      DO i = 1, number_of_right_hand_sides
         dvr_mat(0)%Exact_Solution(1:n_final,i) =  grid%grid_points(:) *                            &
                                           ( grid%grid_points(:) *                                  &
                                             grid%grid_points(:) - one ) * sixth
      END DO
  ELSE IF(type_inhomo == 'exponential') THEN
!
!         Solution is exp(-x) + x *  ( 1 -exp(-1) ) - 1  
!
      added = exp(-one)
      DO i = 1, number_of_right_hand_sides
         dvr_mat(0)%Exact_Solution(1:n_final,i) =  exp(-grid%grid_points(:))                        &
                                                     +                                              &
                                             grid%grid_points(:)                                    &
                                                     *                                              &
                                             ( one - added ) - one                    
      END DO
  ELSE IF(type_inhomo == 'coulomb') THEN      
!
!         Solution is x * ( ln (x) - ln(1) )
!
      added = log(one)
      DO i = 1, number_of_right_hand_sides
         dvr_mat(0)%Exact_Solution(1:n_final,i) =  grid%grid_points(:)                              &
                                                  *                                                 &
                                            ( log ( grid%grid_points(:) ) - added )
      END DO
  ELSE
      Call lnkerr('not valid')
  END IF
!
!     This converts to a solution which has the BC that it is one at the last point
!
  DO i= 1, n_final
     dvr_mat(0)%Exact_Solution(i,:)  = dvr_mat(0)%Exact_Solution(i,:) + grid%grid_points(i)       
  END DO
  dvr_mat(0)%Exact_Solution(n_total,:) = one             
END SUBROUTINE Exact_Cartesian
!***********************************************************************
!***********************************************************************
!deck Exact_Spherical
!***begin prologue     Exact_Spherical
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Exact solutions to the Poisson
!***                   equation for some specific cases.
!***
!***references
!***routines called
!***end prologue       Exact_Spherical
  SUBROUTINE Exact_Spherical(grid, name_spherical)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE(spherical)                          :: name_spherical
  REAL(idp)                                :: scale  
  INTEGER                                  :: info
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  CHARACTER(LEN=3)                         :: itoc
  CHARACTER(LEN=80)                        :: title
!
!         Some tests for cases where solution is one at r = R.  The boundary condition
!         at r=0 is that the solution is a constant.  The constant is specified only
!         by looking at the solution of the inhomogeneous equation as given by the 
!         Green's function.  
!
!         Calculate the solution which is zero at the boundary.
!
  IF (type_inhomo == 'one') THEN
!
!            Solution is (r * r - R * R) / 6 
!
     DO i = 1, number_of_right_hand_sides
        dvr_mat(0)%Exact_solution(1:n_final,i) =  ( grid%grid_points(1:n_final) &
                                                           *                    &
                                                    grid%grid_points(1:n_final) &
                                                           -                    &
                                                         R_max                  &
                                                           *                    &
                                                         R_max )                &
                                                           *                    &
                                                         sixth 
     END DO
!
  ELSE IF(type_inhomo == 'linear') THEN
!
!            Solution is (r * r * r - R * R * R) / 12 
!
     DO i = 1, number_of_right_hand_sides                                     
        dvr_mat(0)%Exact_solution(1:n_final,i) =  ( grid%grid_points(1:n_final) &
                                                       *                        &
                                              grid%grid_points(1:n_final)       &
                                                       *                        &
                                              grid%grid_points(1:n_final)       &
                                                       -                        &
                                                       R_max                    &
                                                        *                       &
                                                       R_max                    &
                                                        *                       &
                                                       R_max ) * sixth * half
     END DO
  ELSE IF(type_inhomo == 'coulomb') THEN      
!
!            Solution is (r - R) / 2 
!
     DO i = 1, number_of_right_hand_sides
        dvr_mat(0)%Exact_solution(1:n_final,i) =  ( grid%grid_points(1:n_final) - R_max )     &
                                                                  *                           &
                                                                half
     END DO
  ELSE
!
     Call lnkerr('not valid')
!
  END IF
!
  DO i = 1, n_final
     dvr_mat(0)%Exact_solution(i,:)    = dvr_mat(0)%Exact_solution(i,:) + R_max
  END DO
  dvr_mat(0)%Exact_solution(n_total,:)    = R_max             
END SUBROUTINE Exact_Spherical
!***********************************************************************
!***********************************************************************
!deck Exact_Spheroidal
!***begin prologue     Exact_Spheroidal
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Exact solutions to the Poisson
!***                   equation for some specific cases.
!***
!***references
!***routines called
!***end prologue       Exact_Spheroidal
  SUBROUTINE Exact_Spheroidal(grid, name_spheroidal)
  IMPLICIT NONE
  TYPE(coordinates)                        :: grid  
  TYPE(spheroidal)                         :: name_spheroidal
  REAL(idp)                                :: scale  
  INTEGER                                  :: info
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  CHARACTER(LEN=3)                         :: itoc
  CHARACTER(LEN=80)                        :: title
  DO i = 1, number_of_right_hand_sides
     dvr_mat(0)%RHS(1:n_final,i) = dvr_mat(0)%RHS(1:n_final,i)          &
                               /                        &
                           sqrt( grid%grid_weights(1:n_final) )
  END DO
END SUBROUTINE Exact_Spheroidal
!***********************************************************************
!***********************************************************************
           END MODULE Poisson_Module
!***********************************************************************
!***********************************************************************
