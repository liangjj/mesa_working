!***********************************************************************
! Grid_Functions
!**begin prologue     Grid_Functions
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           FEDVR, grid, polynomials
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Calculate the regional FEDVR points, weights,
!***                  polynomials and normalizations factors as well as
!***                  some global poibnts and weights.
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***end prologue      Grid_Functions
!***********************************************************************
!***********************************************************************
                           MODULE Grid_Functions
                           USE Renormalization
!***********************************************************************
!***********************************************************************
                              CONTAINS
!**********************************************************************   
!**********************************************************************   
!deck Shape_Functions
!***begin prologue     Shape_Functions
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid, weights and coordinate factors
!***                   needed to construct the KE matrix elements.
!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Shape_Functions(reg)
  IMPLICIT NONE         
  TYPE(regional), DIMENSION(:)                :: reg
  INTEGER                                     :: i
  CHARACTER (LEN=3)                           :: itoc
!
!      Get the points and weights and store them in the q_temp and wt_tmp arrays
!
  DO i = 1, nreg
     ALLOCATE( reg(i)%q(1:reg(i)%n_pts), reg(i)%wt(1:reg(i)%n_pts) )
     Call gauss(reg(i)%q, reg(i)%wt, reg(i)%edge,              & 
                type_quadrature=reg(i)%type_quadrature,        &
                fixed_point=reg(i)%fixed_point, n=reg(i)%n_pts)
  END DO
  IF (prnt(1) == .true.) THEN
      Write(iout,1)
      DO i = 1, nreg
         call Print_Matrix(type_real_vector,reg(i)%q,title='Points Region-'//itoc(i))
         call Print_Matrix(type_real_vector,reg(i)%wt,title='Weights Region-'//itoc(i))
      END DO
  END IF
!
  DO i = 1, nreg 
     ALLOCATE( reg(i)%p(1:reg(i)%n_pts,1:reg(i)%n_pts), reg(i)%dp(1:reg(i)%n_pts,1:reg(i)%n_pts),  &
               reg(i)%ddp(1:reg(i)%n_pts,1:reg(i)%n_pts), reg(i)%q_fac(1:reg(i)%n_pts),            &
               reg(i)%inv_q_fac(1:reg(i)%n_pts), reg(i)%inv_sqrt_q_fac(1:reg(i)%n_pts),            &
               reg(i)%normalization(1:reg(i)%n_pts) )
!
     CALL cpoly(reg(i)%p,reg(i)%dp,reg(i)%ddp,reg(i)%q,reg(i)%n_pts)
!
  END DO
!
  IF (prnt(3) == .true.) THEN
      Write(iout,3)
      DO i = 1, nreg
         call Print_Matrix(type_real_matrix,reg(i)%p,reg(i)%n_pts,reg(i)%n_pts,                     &
                           title='Unnormalized Polynomials Region-'//itoc(i))
         call Print_Matrix(type_real_matrix,reg(i)%dp,reg(i)%n_pts,reg(i)%n_pts,                    &
                           title='First Derivative of Unnormalized Polynomials Region-'//itoc(i))
         call Print_Matrix(type_real_matrix,reg(i)%ddp,reg(i)%n_pts,reg(i)%n_pts,                   &
                           title='Second Derivative of Unnormalized Polynomials Region-'//itoc(i))
     END DO
  END IF  
!
  IF (coordinate_system == 'cartesian') THEN
      Call Coordinate_Functions(grid%xyz,reg)
  ELSE IF(coordinate_system == 'spherical') THEN
      Call Coordinate_Functions(grid%r_theta,reg)
  ELSE IF(coordinate_system == 'cylindrical') THEN
      Call Coordinate_Functions(grid%rho_z,reg)
  ELSE IF(coordinate_system == 'spheroidal') THEN
      Call Coordinate_Functions(grid%xi_eta,reg)
  END IF 
!
!
  IF (prnt(2) == .true.) THEN
      Write(iout,2)
      DO i = 1, nreg
         call Print_Matrix(type_real_vector,reg(i)%q_fac,title='q_factor Region-'//itoc(i))
         call Print_Matrix(type_real_vector,reg(i)%inv_q_fac,title='inverse_q_factor Region-'//itoc(i))
         call Print_Matrix(type_real_vector,reg(i)%inv_sqrt_q_fac,title='inverse_sqrt_q_factor Region-'//itoc(i))
     END DO
  END IF  
  Call Normalization(reg)
  Call Normalized_Functions(reg)
  Call Global_Points_Weights(grid,reg)
!
! Find the first and last function in each region after enforcing boundary conditions
!
  reg(1:nreg)%first = int_one
  reg(1:nreg)%last  = reg(1:nreg)%n_pts
  reg(1:nreg)%n_fun = reg(1:nreg)%n_pts
!
!
  IF (nreg == int_one) THEN
      i = int_one
      IF (grid%drop(1) == .true. ) THEN
          reg(i)%first = int_two
          reg(i)%n_fun = reg(i)%n_fun - int_one
      END IF
      IF (grid%drop(2) == .true. ) THEN
          reg(i)%n_fun = reg(i)%n_fun - int_one
          reg(i)%last = reg(i)%n_fun
      END IF
  ELSE
      i = int_one
      IF (grid%drop(1) == .true. ) THEN
          reg(i)%first = int_two
          reg(i)%n_fun = reg(i)%n_fun - int_one
      END IF
      i = nreg
      IF (grid%drop(2) == .true. ) THEN
          reg(i)%n_fun = reg(i)%n_fun - int_one
          reg(i)%last = reg(i)%n_fun - int_one
      END IF
  END IF
!
!
  IF (prnt(3) == .true.) THEN
      DO i = 1, nreg
         call Print_Matrix(type_real_matrix,reg(i)%bf(:,reg(i)%first:reg(i)%last),reg(i)%n_pts,reg(i)%n_fun,   &
                           title='Final Normalized Basis Polynomials Region-'//itoc(i))
         call Print_Matrix(type_real_matrix,reg(i)%dbf(:,reg(i)%first:reg(i)%last),reg(i)%n_pts,reg(i)%n_fun,  &
                           title='First Derivative of Final Normalized Basis Polynomials Region-'//itoc(i))
         call Print_Matrix(type_real_matrix,reg(i)%ddbf(:,reg(i)%first:reg(i)%last),reg(i)%n_pts,reg(i)%n_fun, &
                           title='Second Derivative of Final Normalized Basis Polynomials Region-'//itoc(i))
     END DO
  END IF  
!
1    FORMAT(/,1x,'DVR Points and Weights') 
2    FORMAT(/,1x,'DVR Coordinate Factors') 
3    FORMAT(/,1x,'Functions and Derivatives before Imposing Boundary Conditions to Remove Functions') 
4    FORMAT(/,1x,'Functions and Derivatives after Removing Functions')
5    Format(/,1x,'Number of Functions in Region = ',i3,' are ',i4)
END SUBROUTINE Shape_Functions
!***********************************************************************
!***********************************************************************
          END MODULE Grid_Functions
!***********************************************************************
!***********************************************************************
