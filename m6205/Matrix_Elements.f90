!***********************************************************************
! Matrix_Elements
!**begin prologue     Matrix_Elements
!**date written       090119   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Iterative, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Calculate the regional kinetic energy matrix elements
!***                  in a FEDVR basis.  
!***description       The routines in this module compute the regional kinetic
!***                  energy matrix elements.  For the FEDVR basis sets 
!***                  based on Gauss quadatures there are cases where one needs 
!***                  to explicitlyxtract a square root singularity before using the FEDVR.
!***                  To do that we define an even and odd type
!***references
!***modules needed    See USE statements below
!***comments          
!***                  
!***                  
!***                  
!***                  
!***end prologue      Matrix_Elements
!***********************************************************************
!***********************************************************************
                           MODULE Matrix_Elements
                           USE Data_Module
                           USE FEDVR_Lib_Shared
                           USE Derived_Types
                           USE FEDVR_Potential
                           USE Renormalization
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!**********************************************************************
!
                            INTERFACE Coordinate_Matrices 
                       MODULE PROCEDURE Radial_Matrices,     &
                                        Theta_Matrices,      &                                             
                                        Phi_Matrices                                                                
                            END INTERFACE Coordinate_Matrices
!
                            INTERFACE Type_Matrices                                                                    
                       MODULE PROCEDURE Even_Matrices,       &
                                        Odd_Matrices                                             
                            END INTERFACE Type_Matrices

!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Matrices
!***begin prologue     Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Matrices(grid, reg, wa )
  IMPLICIT NONE         
  TYPE(coordinates)                                :: grid
  TYPE(regional), DIMENSION(:)                     :: reg
  TYPE (working_arrays)                            :: wa
  TYPE(final_matrices), DIMENSION(:), ALLOCATABLE  :: mat
  TYPE(even_matrices)                              :: even
  TYPE(odd_matrices)                               :: odd
  INTEGER                                          :: i
  CHARACTER(LEN=3)                                 :: itoc
  CHARACTER(LEN=8)                                 :: coordinate_type
!
!
! Set the potential matrix elements
!
  Call PE_DVR_Matrix(reg)
  IF (coordinate_system == 'cartesian') THEN
!
!     Compute either the even, odd or both matrix elements depending on the case
!
      Call Matrix_Elements(even=even, reg=reg)
!     
!     Normalize the matrix
!
      Call Normalized_Matrix_Elements(even, reg)
      lower = int_zero
      upper = int_zero
      skip = int_one
      coordinate_type = grid%xyz%axis
      DO i = 1,nreg
         ALLOCATE( reg(i)%mat(lower:upper) )
      END DO
!
!     Compute the final matrix elements for all cases.
!
     Call Final_matrix_Elements(even=even, reg=reg, mat=mat, type=coordinate_type)
!
!------------------------------------------------------------------------------------
!                 For all the other coordinate systems and components do as above
!                 for cartesian system
!-----------------------------------------------------------------------------------
  ELSE IF(coordinate_system == 'spherical') THEN
      coordinate_type = grid%r_theta%axis
      IF (grid%r_theta%axis == 'r') THEN
          Call Matrix_Elements(even, reg)
          Call Normalized_Matrix_Elements(even, reg)
          lower = int_zero
          upper = grid%r_theta%l_max
          skip = int_one
          DO i = 1,nreg
             ALLOCATE( reg(i)%mat(lower:upper) )
          END DO
          Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
          IF ( grid%drop(1) == .false.) THEN
               Call Transform_Matrix_to_Standard_Form ( reg, mat )
          END IF
      ELSE IF (grid%r_theta%axis == 'theta') THEN
          lower = int_zero
          upper = grid%r_theta%m_max
          skip = int_two
          DO i = 1,nreg
             ALLOCATE( reg(i)%mat(lower:upper) )
          END DO
          Call Matrix_Elements(even, reg)
          Call Normalized_Matrix_Elements(even, reg)
          Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
          IF (grid%r_theta%m_max > 0 ) THEN
              pre_factor = -one
              Call Matrix_Elements(odd,reg)
              Call Normalized_Matrix_Elements(odd,reg)
              lower = int_one
              Call Final_matrix_Elements( odd=odd, reg=reg, mat=mat, type=coordinate_type )
          END IF
          lower = int_zero
          skip = int_one
      END IF
  ELSE IF(coordinate_system == 'cylindrical') THEN
      coordinate_type=grid%rho_z%axis
      Call Matrix_Elements(even, reg)
      Call Normalized_Matrix_Elements(even, reg)
      lower = int_zero
      upper = grid%rho_z%m_max
      skip = int_one
      DO i = 1,nreg
         ALLOCATE( reg(i)%mat(lower:upper) )
      END DO
      Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
      Call Transform_Matrix_to_Standard_Form ( reg, mat )
  ELSE IF(coordinate_system == 'spheroidal') THEN
      coordinate_type=grid%xi_eta%axis
      lower = int_zero
      upper = grid%xi_eta%m_max
      skip = int_two
      IF(grid%xi_eta%axis == 'eta') THEN
         pre_factor = -one
      ELSE IF(grid%xi_eta%axis == 'xi') THEN
         pre_factor = one
      END IF
      Call Matrix_Elements(even, reg)
      Call Normalized_Matrix_Elements(even, reg)
      DO i = 1,nreg
         ALLOCATE( reg(i)%mat(lower:upper) )
      END DO
      Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
      IF (grid%xi_eta%m_max > 0 ) THEN
          Call Matrix_Elements(odd, reg)
          Call Normalized_Matrix_Elements(odd,reg)
          lower = int_one
          Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
          lower = int_zero
          skip = int_one
      END IF
      Call Transform_Matrix_to_Standard_Form ( reg, mat )
  END IF
!
  ALLOCATE( wa%mat(lower:upper) )  
  write(iout,*)
  write(iout,*) '                              Forming the Final Full Matrix'
  write(iout,*)
  Call Form_Global_Matrix(reg, wa, mat )
  IF (prnt(8) == .true. ) THEN
      DO lm = lower, upper, skip
         Call Print_Matrix(type_real_matrix,wa%mat(lm)%ham,nphy,nphy,title='Full Hamiltonian Matrix lm-'//itoc(lm))
      END DO
  END IF
!
  Call Diagonalize_Global_Matrix( wa, mat )   
END SUBROUTINE Matrices
!***********************************************************************
!***********************************************************************
!deck Radial_Matrices
!***begin prologue     Radial_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Radial_Matrices

  SUBROUTINE Radial_Matrices(atom, rad, shl, wa )
  IMPLICIT NONE         
  TYPE(ATOMS)                                      :: atom
  TYPE(RADIAL)                                     :: rad
  TYPE(REGIONAL), DIMENSION(:)                     :: shl
  TYPE(EVEN_MATRICES)                              :: even
  TYPE(FINAL_MATRICES), DIMENSION(:), ALLOCATABLE  :: mat
  TYPE (WORKING_ARRAYS)                            :: wa
  INTEGER                                          :: i
  CHARACTER(LEN=3)                                 :: itoc
!
  Call Type_Matrices(atom, shl, even, atom%n_shl)
  Call Normalized_Matrices(shl, even, atom%n_shl)
  lower = int_zero
  upper = ltop
  skip = int_one
  DO i = 1,atom%n_shl
     ALLOCATE( reg(i)%mat(lower:upper) )
  END DO
  Call Final_matrix_Elements( even=even, reg=reg, mat=mat, n_reg=atom%n_shl, type='r' )
  IF ( grid%drop(1) == .false.) THEN
       Call Transform_Matrix_to_Standard_Form ( reg, mat )
  END IF
!
  ALLOCATE( wa%mat(lower:upper) )  
  write(iout,*)
  write(iout,*) '                              Forming the Final Full Matrix'
  write(iout,*)
  Call Form_Global_Matrix(reg, wa, mat, atom% )
  IF (prnt(3) == .true. ) THEN
      DO lm = lower, upper, skip
         Call Print_Matrix(type_real_matrix,wa%mat(lm)%ham,nphy,nphy,title='Full Hamiltonian Matrix lm-'//itoc(lm))
      END DO
  END IF
!
  Call Diagnalize_Global_Matrix( wa, mat )   
END SUBROUTINE Radial_Matrices
!**********************************************************************
!***********************************************************************
!deck Theta_Matrices
!***begin prologue     Theta_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Theta_Matrices(grid, reg, wa )
  IMPLICIT NONE         
  TYPE(ATOMS)                                      :: atom
  TYPE(THETA)                                      :: theta
  TYPE(REGIONAL), DIMENSION(:)                     :: reg
  TYPE(EVEN_MATRICES)                              :: even
  TYPE(FINAL_MATRICES), DIMENSION(:), ALLOCATABLE  :: mat
  TYPE (WORKING_ARRAYS)                            :: wa
  INTEGER                                          :: i
  CHARACTER(LEN=3)                                 :: itoc
  CHARACTER(LEN=8)                                 :: coordinate_type
!
!
! Set the potential matrix elements
!
  Call PE_DVR_Matrix(reg)
  IF (coordinate_system == 'cartesian') THEN
!
!     Compute either the even, odd or both matrix elements depending on the case
!
      Call Matrix_Elements(even=even, reg=reg)
!     
!     Normalize the matrix
!
      Call Normalized_Matrix_Elements(even, reg)
      lower = int_zero
      upper = int_zero
      skip = int_one
      coordinate_type = grid%xyz%axis
      DO i = 1,nreg
         ALLOCATE( reg(i)%mat(lower:upper) )
      END DO
!
!     Compute the final matrix elements for all cases.
!
     Call Final_matrix_Elements(even=even, reg=reg, mat=mat, type=coordinate_type)
!
!------------------------------------------------------------------------------------
!                 For all the other coordinate systems and components do as above
!                 for cartesian system
!-----------------------------------------------------------------------------------
  ELSE IF(coordinate_system == 'spherical') THEN
      coordinate_type = grid%r_theta%axis
      IF (grid%r_theta%axis == 'r') THEN
          Call Matrix_Elements(even, reg)
          Call Normalized_Matrix_Elements(even, reg)
          lower = int_zero
          upper = grid%r_theta%l_max
          skip = int_one
          DO i = 1,nreg
             ALLOCATE( reg(i)%mat(lower:upper) )
          END DO
          Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
          IF ( grid%drop(1) == .false.) THEN
               Call Transform_Matrix_to_Standard_Form ( reg, mat )
          END IF
      ELSE IF (grid%r_theta%axis == 'theta') THEN
          lower = int_zero
          upper = grid%r_theta%m_max
          skip = int_two
          DO i = 1,nreg
             ALLOCATE( reg(i)%mat(lower:upper) )
          END DO
          Call Matrix_Elements(even, reg)
          Call Normalized_Matrix_Elements(even, reg)
          Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
          IF (grid%r_theta%m_max > 0 ) THEN
              pre_factor = -one
              Call Matrix_Elements(odd,reg)
              Call Normalized_Matrix_Elements(odd,reg)
              lower = int_one
              Call Final_matrix_Elements( odd=odd, reg=reg, mat=mat, type=coordinate_type )
          END IF
          lower = int_zero
          skip = int_one
      END IF
  ELSE IF(coordinate_system == 'cylindrical') THEN
      coordinate_type=grid%rho_z%axis
      Call Matrix_Elements(even, reg)
      Call Normalized_Matrix_Elements(even, reg)
      lower = int_zero
      upper = grid%rho_z%m_max
      skip = int_one
      DO i = 1,nreg
         ALLOCATE( reg(i)%mat(lower:upper) )
      END DO
      Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
      Call Transform_Matrix_to_Standard_Form ( reg, mat )
  ELSE IF(coordinate_system == 'spheroidal') THEN
      coordinate_type=grid%xi_eta%axis
      lower = int_zero
      upper = grid%xi_eta%m_max
      skip = int_two
      IF(grid%xi_eta%axis == 'eta') THEN
         pre_factor = -one
      ELSE IF(grid%xi_eta%axis == 'xi') THEN
         pre_factor = one
      END IF
      Call Matrix_Elements(even, reg)
      Call Normalized_Matrix_Elements(even, reg)
      DO i = 1,nreg
         ALLOCATE( reg(i)%mat(lower:upper) )
      END DO
      Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
      IF (grid%xi_eta%m_max > 0 ) THEN
          Call Matrix_Elements(odd, reg)
          Call Normalized_Matrix_Elements(odd,reg)
          lower = int_one
          Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
          lower = int_zero
          skip = int_one
      END IF
      Call Transform_Matrix_to_Standard_Form ( reg, mat )
  END IF
!
  ALLOCATE( wa%mat(lower:upper) )  
  write(iout,*)
  write(iout,*) '                              Forming the Final Full Matrix'
  write(iout,*)
  Call Form_Global_Matrix(reg, wa, mat )
  IF (prnt(8) == .true. ) THEN
      DO lm = lower, upper, skip
         Call Print_Matrix(type_real_matrix,wa%mat(lm)%ham,nphy,nphy,title='Full Hamiltonian Matrix lm-'//itoc(lm))
      END DO
  END IF
!
  Call Diagonalize_Global_Matrix( wa, mat )   
END SUBROUTINE Theta_Matrices
!***********************************************************************
!***********************************************************************
!deck Matrices
!***begin prologue     Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Matrices(grid, reg, wa )
  IMPLICIT NONE         
  TYPE(coordinates)                                :: grid
  TYPE(regional), DIMENSION(:)                     :: reg
  TYPE (working_arrays)                            :: wa
  TYPE(final_matrices), DIMENSION(:), ALLOCATABLE  :: mat
  TYPE(even_matrices)                              :: even
  TYPE(odd_matrices)                               :: odd
  INTEGER                                          :: i
  CHARACTER(LEN=3)                                 :: itoc
  CHARACTER(LEN=8)                                 :: coordinate_type
!
!
! Set the potential matrix elements
!
  Call PE_DVR_Matrix(reg)
  IF (coordinate_system == 'cartesian') THEN
!
!     Compute either the even, odd or both matrix elements depending on the case
!
      Call Matrix_Elements(even=even, reg=reg)
!     
!     Normalize the matrix
!
      Call Normalized_Matrix_Elements(even, reg)
      lower = int_zero
      upper = int_zero
      skip = int_one
      coordinate_type = grid%xyz%axis
      DO i = 1,nreg
         ALLOCATE( reg(i)%mat(lower:upper) )
      END DO
!
!     Compute the final matrix elements for all cases.
!
     Call Final_matrix_Elements(even=even, reg=reg, mat=mat, type=coordinate_type)
!
!------------------------------------------------------------------------------------
!                 For all the other coordinate systems and components do as above
!                 for cartesian system
!-----------------------------------------------------------------------------------
  ELSE IF(coordinate_system == 'spherical') THEN
      coordinate_type = grid%r_theta%axis
      IF (grid%r_theta%axis == 'r') THEN
          Call Matrix_Elements(even, reg)
          Call Normalized_Matrix_Elements(even, reg)
          lower = int_zero
          upper = grid%r_theta%l_max
          skip = int_one
          DO i = 1,nreg
             ALLOCATE( reg(i)%mat(lower:upper) )
          END DO
          Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
          IF ( grid%drop(1) == .false.) THEN
               Call Transform_Matrix_to_Standard_Form ( reg, mat )
          END IF
      ELSE IF (grid%r_theta%axis == 'theta') THEN
          lower = int_zero
          upper = grid%r_theta%m_max
          skip = int_two
          DO i = 1,nreg
             ALLOCATE( reg(i)%mat(lower:upper) )
          END DO
          Call Matrix_Elements(even, reg)
          Call Normalized_Matrix_Elements(even, reg)
          Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
          IF (grid%r_theta%m_max > 0 ) THEN
              pre_factor = -one
              Call Matrix_Elements(odd,reg)
              Call Normalized_Matrix_Elements(odd,reg)
              lower = int_one
              Call Final_matrix_Elements( odd=odd, reg=reg, mat=mat, type=coordinate_type )
          END IF
          lower = int_zero
          skip = int_one
      END IF
  ELSE IF(coordinate_system == 'cylindrical') THEN
      coordinate_type=grid%rho_z%axis
      Call Matrix_Elements(even, reg)
      Call Normalized_Matrix_Elements(even, reg)
      lower = int_zero
      upper = grid%rho_z%m_max
      skip = int_one
      DO i = 1,nreg
         ALLOCATE( reg(i)%mat(lower:upper) )
      END DO
      Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
      Call Transform_Matrix_to_Standard_Form ( reg, mat )
  ELSE IF(coordinate_system == 'spheroidal') THEN
      coordinate_type=grid%xi_eta%axis
      lower = int_zero
      upper = grid%xi_eta%m_max
      skip = int_two
      IF(grid%xi_eta%axis == 'eta') THEN
         pre_factor = -one
      ELSE IF(grid%xi_eta%axis == 'xi') THEN
         pre_factor = one
      END IF
      Call Matrix_Elements(even, reg)
      Call Normalized_Matrix_Elements(even, reg)
      DO i = 1,nreg
         ALLOCATE( reg(i)%mat(lower:upper) )
      END DO
      Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
      IF (grid%xi_eta%m_max > 0 ) THEN
          Call Matrix_Elements(odd, reg)
          Call Normalized_Matrix_Elements(odd,reg)
          lower = int_one
          Call Final_matrix_Elements( even=even, reg=reg, mat=mat, type=coordinate_type )
          lower = int_zero
          skip = int_one
      END IF
      Call Transform_Matrix_to_Standard_Form ( reg, mat )
  END IF
!
  ALLOCATE( wa%mat(lower:upper) )  
  write(iout,*)
  write(iout,*) '                              Forming the Final Full Matrix'
  write(iout,*)
  Call Form_Global_Matrix(reg, wa, mat )
  IF (prnt(8) == .true. ) THEN
      DO lm = lower, upper, skip
         Call Print_Matrix(type_real_matrix,wa%mat(lm)%ham,nphy,nphy,title='Full Hamiltonian Matrix lm-'//itoc(lm))
      END DO
  END IF
!
  Call Diagonalize_Global_Matrix( wa, mat )   
END SUBROUTINE Matrices
!***********************************************************************
!***********************************************************************
!deck Even_Matrices
!***begin prologue     KE
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate the parts of the regional matrices
!***                   not depending on the (l,m) quantum numbers.
!***                   Also, the formulas programmed assume that 1) we are dealing with Del^2
!***                   left in its standard form and no first derivative have been transformed
!***                   away and 2) the volume element is included in the definition of the the matrix
!***                   element.  There is one special case.  The radial problem in spherical coordinates
!***                   may be treated using the cartesian form but requiring the functions to vanish at the
!***                   origin.  This was done for consistency with older code. The programmed expression 
!***                   is what is obtained after integrating by parts.  if there is a singular point at 
!***                   the boundaries, the integrand remains well behaved and it is not necessary to 
!***                   make that a quadrature point.  This allows us to use Gauss-Radau rules for those elements instead
!***                   of Gauss-Lobatto and you do not have to remove the first basis function.  This was
!***                   first pointed out to me by Brett Esry.  This routine will handle cases where
!***                   the FEDVR basis functions are polynomials.  Finally, the integrand has been integrated by parts, 
!***                   which results in the minus sign below.
!***references

!***routines called    iosys, util and mdutil
!***end prologue       Even_Matrices
  SUBROUTINE Even_Matrices(atom,reg,even,n_reg)
  IMPLICIT NONE
  TYPE(REGIONAL), DIMENSION(:)         :: reg
  TYPE(EVEN_MATRICES)                  :: even
  INTEGER                              :: n_reg
  INTEGER                              :: ir
  INTEGER                              :: i
  INTEGER                              :: k
  CHARACTER(LEN=3)                     :: itoc
!
!
!
  DO i = 1, n_reg
     ALLOCATE( reg(i)%even%tr(1:reg(i)%n_pts,1:reg(i)%n_pts) )
     reg(i)%even%tr(:,:) = zero
  END DO
!
!    Loop over regions
!
  DO ir = 1, n_reg
!
!    Loop over (i,j) 
!
!
     DO i = 1, reg(ir)%n_pts
!
!       Sum over quadrature points
!
        DO k = 1, reg(ir)%n_pts
           reg(ir)%even%tr(i,1:i) = reg(ir)%even%tr(i,1:i)   &
                                  -                          &
                 reg(ir)%q_fac(k) * reg(ir)%wt(k)            &
                                  *                          &
                 reg(ir)%dp(k,i)  * reg(ir)%dp(k,1:i) 
        END DO
        reg(ir)%even%tr(1:i,i) = reg(ir)%even%tr(i,1:i)
     END DO
  END DO
  IF (prnt(2) == .true. ) THEN
      DO i = 1, n_reg
         call Print_Matrix(type_real_matrix,reg(i)%even%tr,reg(i)%n_pts,reg(i)%n_pts,     &
                           title='Unnormalized_Even_Nabla_Matrix-Region-'//itoc(i))
      END DO
  END IF
!
END SUBROUTINE Even_Matrices
!***********************************************************************
!***********************************************************************
!deck Odd_Matrices
!***begin prologue     Odd_Matrices
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            This routine will handle FEDVR basis functions when there
!***                   is a square root singularity that needs to be removed explicitly.
!***                   Once the square root factor is removed the remaining part of
!***                   each basis function may be taken to be of polynomial form.
!***                   This is needed for certain variables such as the odd m legendre
!***                   or spheroidal functions.  Again, the programmed formulas are
!***                   obtained after integration by parts.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Odd_KE

  SUBROUTINE Odd_Matrices(atom,reg,odd,n_reg)
  IMPLICIT NONE
  TYPE(REGIONAL), DIMENSION(:)         :: reg
  TYPE(ODD_MATRICES)                   :: odd
  INTEGER                              :: n_reg
  INTEGER                              :: ir
  INTEGER                              :: i
  INTEGER                              :: k
  CHARACTER(LEN=3)                     :: itoc
!
!
  DO i = 1, n_reg
     ALLOCATE( reg(i)%odd%tr(1:reg(i)%n_pts,1:reg(i)%n_pts) )
     reg(i)%odd%tr(:,:) = zero
  END DO
!
!    Loop over the regions
!
  DO ir = 1, n_reg
!
!    Loop over the (i,j) with j an implicit loop
!

     DO  i = 1, reg(ir)%n_pts
!
!        Sum over the quadrature points
!
         DO k = 1, reg(ir)%n_pts
            reg(ir)%odd%tr(i,1:i) = reg(ir)%odd%tr(i,1:i)    &
                                  -                          &
                 reg(ir)%q_fac(k) * reg(ir)%q_fac(k)         &
                                  *                          &
                 reg(ir)%wt(k)    * reg(ir)%dp(k,i)          &
                                  * reg(ir)%dp(k,1:i)
         END DO
         reg(ir)%odd%tr(i,1:i) = reg(ir)%odd%tr(i,1:i)       &
                               -                             &
                    pre_factor                               &
                               *                             &
                   ( reg(ir)%q(i)                            &
                               *                             &
                     reg(ir)%q_fac(i)                        &
                               *                             &
                     reg(ir)%p(i,i)                          &
                               *                             &
                     reg(ir)%dp(i,1:i) )                     &
                               *                             &
                     reg(ir)%wt(i)
         DO k = 1 , i
            reg(ir)%odd%tr(i,k) = reg(ir)%odd%tr(i,k)        &
                                -                            &
                     pre_factor                              &
                                *                            &
                     ( reg(ir)%q(k)                          &
                                *                            &
                       reg(ir)%q_fac(k)                      &
                                *                            &
                       reg(ir)%p(k,k)                        &
                                *                            &
                       reg(ir)%dp(k,i) )                     &
                                *                            &
                       reg(ir)%wt(k)         
         END DO
         reg(ir)%odd%tr(i,i)  = reg(ir)%odd%tr(i,i)          &
                              -                              &
         reg(ir)%q(i)         * reg(ir)%q(i)                 &
                              *                              &
         reg(ir)%p(i,i)       * reg(ir)%p(i,i)               &
                              *                              &
                       reg(ir)%wt(i)         
     END DO
!
!    Symmetrize
!
     DO  i = 1, reg(ir)%n_pts
         reg(ir)%odd%tr(i,1:i) =                             &
                reg(ir)%inv_sqrt_q_fac(i)                    &
                           *                                 &
                reg(ir)%odd%tr(i,1:i)                        &
                           *                                 &
                reg(ir)%inv_sqrt_q_fac(1:i)
         reg(ir)%odd%tr(1:i,i) = reg(ir)%odd%tr(i,1:i)
     END DO
  END DO
!
  IF (prnt(2) == .true. ) THEN
      DO i = 1, n_reg
         call Print_Matrix(type_real_matrix,reg(i)%odd%tr,reg(i)%n_pts,reg(i)%n_pts,                           &
                           title='Unnormalized_Odd_Nabla_Matrix-Region-'//itoc(i))
      END DO
  END IF
!
END SUBROUTINE Odd_Matrices
!***********************************************************************
!***********************************************************************
!deck Final_Matrix_Elements
!***begin prologue     Final_Matrix_Elements
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Add angular momentum terms and any potentials
!***references

!***routines called    iosys, util and mdutil
!***end prologue       Final_Matrix_Elements
  SUBROUTINE Final_Matrix_Elements(even, odd, reg, mat, type, n_reg)
  IMPLICIT NONE
  TYPE(REGIONAL), DIMENSION(:)                     :: reg
  TYPE(EVEN_MATRICES), OPTIONAL                    :: even
  TYPE(ODD_MATRICES), OPTIONAL                     :: odd
  TYPE(FINAL_MATRICES), DIMENSION(:)               :: mat
  INTEGER                                          :: n_reg
  INTEGER                                          :: i
  INTEGER                                          :: lm
  CHARACTER(LEN=*)                                 :: type
  CHARACTER(LEN=3)                                 :: itoc
!
!
  dscale= - hbar * hbar * half / mass  
  IF ( PRESENT(even) == .true. ) THEN
       DO i = 1, n_reg
          reg(i)%even%ham(:,:) = dscale * reg(i)%even%ham(:,:) 
          Call Add_Potential( reg(i)%even%ham, reg(i)%pot, reg(i)%q_fac, reg(i)%n_pts )
          DO lm = lower, upper, skip
             ALLOCATE( reg(i)%mat(lm)%ham(1:reg(i)%n_pts,1:reg(i)%n_pts) )
             reg(i)%mat(lm)%ham(:,:) = reg(i)%even%ham(:,:) 
             Call Add_Angular_Momentum ( reg(i)%mat(lm)%ham, reg(i)%q, reg(i)%q_fac, reg(i)%inv_q_fac,  &
                                         lm, type, reg(i)%n_pts )
          END DO
       END DO
  END IF 
  IF ( PRESENT(odd) == .true. ) THEN
       DO i = 1, n_reg
          reg(i)%odd%ham(:,:) = dscale * reg(i)%odd%ham(:,:) 
          Call Add_Potential( reg(i)%odd%ham, reg(i)%pot, reg(i)%q_fac, reg(i)%n_pts )
          DO lm = lower, upper, skip
             ALLOCATE( reg(i)%mat(lm)%ham(1:reg(i)%n_pts,1:reg(i)%n_pts) )
             reg(i)%mat(lm)%ham(:,:) = reg(i)%odd%ham(:,:) 
             Call Add_Angular_Momentum ( reg(i)%mat(lm)%ham, reg(i)%q, reg(i)%q_fac, reg(i)%inv_q_fac,  &
                                         lm, type, reg(i)%n_pts )
          END DO
       END DO
  END IF 
  IF (prnt(4) == .true. ) THEN
      DO i = 1, n_reg
         DO lm = lower, upper, skip
            Call Print_Matrix(type_real_matrix,reg(i)%mat(lm)%ham(reg(i)%first:reg(i)%last,reg(i)%first:reg(i)%last),     &
                              reg(i)%n_fun,reg(i)%n_fun,                                                                  &
                              title='Normalized Hamiltonian Matrix lm-'//itoc(lm)//' Region-'//itoc(i))
         END DO
      END DO
  END IF
END SUBROUTINE Final_Matrix_Elements
!***********************************************************************
!***********************************************************************
!deck Add_Angular_Momentum
!***begin prologue     Add_Angular_Momentum
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Add in the angular momentum.  Note the minus sign.  This is consistent
!***                   with the definition of Del^2 and the integration by parts.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Add_Angular_Momentum
  SUBROUTINE Add_Angular_Momentum( mat, q, q_fac, inv_q_fac, lm, type, n )
  IMPLICIT NONE
  REAL(idp), DIMENSION(:,:)           :: mat
  REAL(idp), DIMENSION(:)             :: q
  REAL(idp), DIMENSION(:)             :: q_fac
  REAL(idp), DIMENSION(:)             :: inv_q_fac
  INTEGER                             :: lm
  CHARACTER (LEN=*)                   :: type
  INTEGER                             :: n
  INTEGER                             :: num
  INTEGER                             :: i
!
  IF ( type == 'r' ) THEN
       num = lm * ( lm + int_one )
       DO i = 1, n
          mat(i,i) = mat(i,i) - dscale *num / ( q(i) * q(i) ) *q_fac(i)
       END DO
  END IF
!
!
  IF ( type == 'theta' .or. type == 'eta') THEN
       num = lm * lm
       DO i = 1, n
          mat(i,i) = mat(i,i) - dscale * num * inv_q_fac(i)
       END DO 
  END IF
END SUBROUTINE Add_Angular_Momentum
!***********************************************************************
!***********************************************************************
!deck Add_Potential
!***begin prologue     Add_Potential
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Add in the potential.  The minus sign is to make it consistent
!***                   with the definition of Del^2 and the integration by parts.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Add_Potential
  SUBROUTINE Add_Potential( mat, pot, q_fac, n )
  IMPLICIT NONE
  REAL(idp), DIMENSION(:,:)                        :: mat
  REAL(idp), DIMENSION(:)                          :: pot
  REAL(idp), DIMENSION(:)                          :: q_fac
  INTEGER                                          :: n
  INTEGER                                          :: i
!
  DO i = 1, n
     mat(i,i) = mat(i,i) + pot(i) * q_fac(i)
  END DO
END SUBROUTINE Add_Potential
!***********************************************************************
!***********************************************************************
!deck Transform_Matrix_to_Standard_Form
!***begin prologue     Transform_Matrix_to_Standard_Form
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            In certain coordinate systems there is a simple metric
!                      due to the structure of the kinetic energy matrix
!                      which turn a standard eigenvalue problem to what looks
!                      like a general eigenvalue problem.  This is trivially
!                      removed by a simple diagonal pre and post multiplication.
!                      This unitary transformation also makes the definition of
!                      the vectors different and this needs to be accounted for later.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE Transform_Matrix_to_Standard_Form(reg, mat )
  IMPLICIT NONE
  TYPE(REGIONAL), DIMENSION(:)                   :: reg
  TYPE(FINAL_MATRICES), DIMENSION(:)             :: mat
  INTEGER                                        :: i
  INTEGER                                        :: j
  CHARACTER(LEN=3)                               :: itoc
!
  Write(iout,*)
  write(iout,*) 'Forming Scaled Hamiltonian'
  DO i = 1, n_reg
     DO lm = lower, upper, skip
        DO j = 1, reg(i)%n_pts
!
!          Pre and post multiplication
!
           reg(i)%mat(lm)%ham(j,1:j) =                                       &
                                       reg(i)%inv_sqrt_q_fac(j)              &
                                               *                             &
                                       reg(i)%mat(lm)%ham(j,1:j)             &
                                               *                             &
                                       reg(i)%inv_sqrt_q_fac(1:j)
           reg(i)%mat(lm)%ham(1:j,j) = reg(i)%mat(lm)%ham(j,1:j) 
        END DO
     END DO
  END DO
!
!
  IF (prnt(4) == .true. ) THEN
      DO i = 1, n_reg
         DO lm = lower, upper, skip
            call Print_Matrix(type_real_matrix,reg(i)%mat(lm)%ham(reg(i)%first:reg(i)%last,reg(i)%first:reg(i)%last),  &
                              reg(i)%n_fun,reg(i)%n_fun,                                                               &
                              title='Scaled Hamiltonian Matrix lm-'//itoc(lm)//' Region-'//itoc(i))
         END DO
      END DO
  END IF
  END SUBROUTINE Transform_Matrix_to_Standard_Form
!***********************************************************************
!***********************************************************************
!deck Form_Global_Matrix.f
!***begin prologue     Form_Global_Matrix
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Form the global matrices.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       
!
  SUBROUTINE Form_Global_Matrix(reg, wa, mat, n_reg )
  IMPLICIT NONE
  TYPE (WORKING_ARRAYS)                   :: wa
  TYPE(REGIONAL), DIMENSION(:)            :: reg
  TYPE(FINAL_MATRICES), DIMENSION(:)      :: mat
  INTEGER                                 :: n_reg
  INTEGER                                 :: i
  INTEGER                                 :: first
  INTEGER                                 :: last
!  
  DO lm = lower, upper, skip
     ALLOCATE ( wa%mat(lm)%ham(1:nphy,1:nphy) )
     wa%mat(lm)%ham(:,:) = 0.d0
     first = 1
     DO i = 1, n_reg
        last = first + reg(i)%n_fun - 1
        wa%mat(lm)%ham(first:last,first:last) = reg(i)%mat(lm)%ham(reg(i)%first:reg(i)%last,reg(i)%first:reg(i)%last)      
        first = last
        DEALLOCATE(reg(i)%mat(lm)%ham)
     END DO    
  END DO
END SUBROUTINE Form_Global_Matrix
!***********************************************************************
!***********************************************************************
!deck Diagonalize_Global_Matrix
!***begin prologue     Diagonalize_Global_Matrix
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Diagonalize the global matrix
!***                   
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Diagonalize_Global_Matrix
!
  SUBROUTINE Diagonalize_Global_Matrix( wa, mat )
  IMPLICIT NONE
  TYPE(FINAL_MATRICES), DIMENSION(:)             :: mat
  TYPE(WORKING_ARRAYS)                           :: wa
  INTEGER                                        :: count
  INTEGER                                        :: i
  INTEGER                                        :: j
  INTEGER                                        :: info
  CHARACTER(LEN=3)                               :: itoc
!  
  write(iout,*)
  write(iout,*) '                              Diagonalizing the Full Matrix'
  write(iout,*)
  ALLOCATE(wa%eigenvectors(1:nphy,1:nphy),                &
           wa%eigenvalues(1:nphy),                        &
           wa%work(3*nphy),                               &
           wa%lower_mat( nphy*(nphy+1)/2 ) )
  DO lm = lower, upper, skip
     count = 0
     DO i = 1, nphy
        DO j = 1, i
           count = count + 1
           wa%lower_mat(count) = wa%mat(lm)%ham(i,j)
        END DO
     END DO
!
!    Used packed form of diagonalizer
!
     Call dspev('v','u',nphy,wa%lower_mat,wa%eigenvalues,wa%eigenvectors,nphy,wa%work,info)
     IF (prnt(10) == .true. ) THEN
         Call Print_Matrix(type_real_vector,wa%eigenvalues,title='Eigenvalues lm-'//itoc(lm))
     END IF
     IF (prnt(11) == .true. ) THEN
         Call Print_Matrix(type_real_matrix,wa%eigenvectors,nphy,nphy,title='Eigenvectors lm-'//itoc(lm))
     END IF
  END DO
  DEALLOCATE(wa%eigenvectors, wa%eigenvalues, wa%work, wa%lower_mat )
END SUBROUTINE Diagonalize_Global_Matrix
!***********************************************************************
           END MODULE Matrix_Elements
!***********************************************************************
