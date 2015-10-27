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
                           USE Data
                           USE Grid_Defined_Types
                           USE Matrix_Print
                           USE Renormalization
!***********************************************************************
!***********************************************************************
!                          Explicit Interfaces
!***********************************************************************
!**********************************************************************
!
                            INTERFACE Coordinate_Matrices 
                       MODULE PROCEDURE Radial_Mat,     &
                                        Theta_Mat,      &                                             
                                        Phi_Mat                                                                
                            END INTERFACE Coordinate_Matrices
!
                            INTERFACE Type_Matrices                                                                    
                       MODULE PROCEDURE Even_Mat,       &
                                        Odd_Mat                                             
                            END INTERFACE Type_Matrices
!
!***********************************************************************
!***********************************************************************
                              CONTAINS
!***********************************************************************
!***********************************************************************
!deck Potential_Mat
!***begin prologue     Potential_Mat
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Calculate the local potential matrix for a number of
!***                   fairly standard potentials.
!***                   
!***                   
!***                   
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       Potential_Mat
  SUBROUTINE Potential_Mat(rad, shl, n_reg)
  IMPLICIT NONE
  TYPE(RADIAL)                     :: rad
  TYPE(REGIONAL), DIMENSION(:)     :: shl
  INTEGER                          :: n_reg
  INTEGER                          :: i
  CHARACTER (LEN=80)               :: type
  REAL(idp)                        :: fpkey
  LOGICAL                          :: dollar
  LOGICAL                          :: logkey
  INTEGER                          :: intkey
  CHARACTER(LEN=80)                :: chrkey
  CHARACTER(LEN=8)                 :: itoc
!
!
  IF ( dollar('$potential',card,cpass,inp) )THEN
       type = chrkey (card,'type_potential','none',' ')
  ELSE
       call lnkerr('Potential Keywords Absent:Quit')
  END If

!
!
  DO i = 1, n_reg
     ALLOCATE( shl(i)%pot(1:shl(i)%n_pts) )
     shl(i)%pot(:) = zero     
!
     IF(type == 'none') then
!
        call none(shl(i)%pot,                                       &
                  shl(i)%q,                                         &
                  shl(i)%n_pts,                                     &
                  prnt)
!
     ELSE IF(type == 'well') then
!
        depth=fpkey(card,'well_depth',zero,' ')
        call vwell(shl(i)%pot,                                      &
                   depth,                                           &
                   shl(i)%n_pts,                                    &
                   prnt)
!
     ELSE IF(type == 'exponential') then
!
        amp(1)=fpkey(card,'amplitude',-1.d0,' ')
        expnt(1)=fpkey(card,'exponent',expnt,' ')
        call vexp(shl(i)%pot,                                       &
                  shl(i)%q,                                         &
                  amp,                                              &
                  expnt,                                            &
                  shl(i)%n_pts,                                     &
                  prnt)
!
     ELSE IF(type == 'yukawa') then
!
        amp(1)=fpkey(card,'amplitude',-1.d0,' ')
        expnt(1)=fpkey(card,'exponent',expnt,' ')
        call vyukawa(shl(i)%pot,                                    &
                     shl(i)%q,                                      &
                     amp,                                           &
                     expnt,                                         &
                     shl(i)%n_pts,                                  &
                     prnt)
!
     ELSE IF(type == 'power_exponential') then
!
        amp(1)=fpkey(card,'amplitude',1.d0,' ')
        expnt(1)=fpkey(card,'exponent',expnt,' ')
        n_p=intkey(card,'power',0,' ')
        call v_pow_exp(shl(i)%pot,                                  &
                       shl(i)%q,                                    &
                       amp,                                         &
                       expnt,                                       &
                       n_p,                                         &
                       shl(i)%n_pts,                                &
                       prnt)
!
     ELSE IF(type == 'sum_exponential') then
!
        call fparr(card,'amplitudes',amp,2,' ')
        call fparr(card,'exponents',expnt,2,' ')
        call vexp_sum(shl(i)%pot,                                   &
                      shl(i)%q,                                     &
                      amp,                                          &
                      expnt,                                        &
                      shl(i)%n_pts,                                 &
                      prnt)
!
     ELSE IF(type == 'coulomb') then
!
        charge=fpkey(card,'charge',-one,' ')
        call vcoul(shl(i)%pot,                                      &
                   shl(i)%q,                                        &
                   charge,                                          &
                   shl(i)%n_pts,                                    &
                   prnt)
!
     ELSE IF(type == 'eberlonium') then
!
        charge=fpkey(card,'charge',-one,' ')
        n_p=intkey(card,'power',0,' ')
        amp(1)=fpkey(card,'a',1.d0,' ')
        amp(2)=fpkey(card,'b',1.d0,' ')
        call v_eberlonium(shl(i)%pot,                               &
                          shl(i)%q,                                 &
                          charge,                                   &
                          amp(1),                                   &
                          amp(2),                                   &
                          n_p,                                      &
                          shl(i)%n_pts,                             &
                          prnt)
!
     ELSE IF(type == 'inverse_r4') then
!
        call vir4(shl(i)%pot,                                       &
                  shl(i)%q,                                         &
                  shl(i)%n_pts,                                     &
                  prnt)
!
     ELSE IF(type == 'rounded_well') then
!
        nwell=intkey(card,'nwell',ten,' ')
        awell=fpkey(card,'a_well',14.d0,' ')
        call vrwell(shl(i)%pot,                                     &
                    shl(i)%q,                                       &
                    awell,                                          &
                    nwell,                                          &
                    shl(i)%n_pts,                                   &
                    prnt)
!
     ELSE IF(type == 'harmonic_oscillator') then
!
!        enter the mass and frequency in atomic units.
!         
        mass=fpkey(card,'mass',1.d0,' ')
        omega=fpkey(card,'omega',1.d0,' ')
        factor= mass * omega * omega * half
        hbar=one
        write(iout,1) mass, omega
        call vhmo(shl(i)%pot,                                       &
                  shl(i)%q,                                         &
                  factor,                                           &
                  shl(i)%n_pts,                                     &
                  prnt)
!
     ELSE IF(type == 'anharmonic_oscillator') then
!
        call vanhmo(shl(i)%pot,                                     &
                    shl(i)%q,                                       &
                    shl(i)%n_pts,                                   &
                    prnt)
!
     ELSE IF(type == 'expres') then
!
        call fparr(card,'amplitude',amp,2,' ')
        call fparr(card,'exponent',expnt,2,' ')
        shift=fpkey(card,'exponent_shift',zero,' ')
        call vres(shl(i)%pot,                                       &
                  shl(i)%q,                                         &
                  amp,                                              &
                  expnt,                                            &
                  shift,                                            &
                  shl(i)%n_pts,                                     &
                  prnt)
!
     ELSE IF(type == 'periodic') then
!
        n_scale=intkey(card,'n_i',10,' ')         
        e_c=fpkey(card,'e_c',.001d0,' ')         
        awell=n_scale/e_c
        call vperiod(shl(i)%pot,                                    &
                     shl(i)%q,                                      &
                     awell,                                         &
                     shl(i)%n_pts,                                  &
                     prnt)
!
     ELSE
        Call lnkerr('screwed up potential. quit')
     END IF
!
  END DO
  IF(prnt(7)) then
     DO i = 1, n_reg
        write(iout,2) type
         call Print_Matrix(type_real_vector,shl(i)%pot,title='potential-shlion-'//itoc(i))
     END DO
  END IF
!
!
 1 Format(/,1x,'oscillator mass      = ',e15.8, &
          /,1x,'oscillator-frequency = ',e15.8)
 2 FORMAT(/,20x,'Potential Type = ',a32)
!
END SUBROUTINE Potential_Mat
!***********************************************************************
!***********************************************************************
!deck Radial_Mat
!***begin prologue     Radial_Mat
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       Radial_Mat

  SUBROUTINE Radial_Mat(atom, rad, shl, wa )
  IMPLICIT NONE         
  TYPE(ATOMS)                                      :: atom
  TYPE(RADIAL)                                     :: rad
  TYPE(REGIONAL), DIMENSION(:)                     :: shl
  TYPE(BOUNDARY_CONDITIONS)                        :: bc
  TYPE(EVEN_MATRICES)                              :: even
  TYPE(FINAL_MATRICES), DIMENSION(:), ALLOCATABLE  :: mat
  TYPE (WORKING_ARRAYS)                            :: wa
  INTEGER                                          :: i
  INTEGER                                          :: lm
  CHARACTER(LEN=3)                                 :: itoc
!
  Call Potential_Matrix(rad, shl, atom%n_shl)
  Call Type_Matrices(shl, even, atom%n_shl)
  Call Normalized_Matrices(shl, even, atom%n_shl)
  lower = int_zero
  upper = ltop
  skip = int_one
  DO i = 1,atom%n_shl
     ALLOCATE( shl(i)%mat(lower:upper) )
  END DO
  Call Final_matrix_Elements( even=even, reg=shl, mat=mat, n_reg=atom%n_shl, type='r' )
  IF ( rad%bc%drop(1) == .false.) THEN
       Call Transform_Matrix_to_Standard_Form ( shl, mat, atom%n_shl )
  END IF
!
  ALLOCATE( wa%mat(lower:upper) )  
  write(iout,*)
  write(iout,*) '                              Forming the Final Full Matrix'
  write(iout,*)
  Call Form_Global_Matrix( shl, wa, mat, atom%n_shl, rad%n_phy )
  IF (prnt(3) == .true. ) THEN
      DO lm = lower, upper, skip
         Call Print_Matrix(type_real_matrix,wa%mat(lm)%ham,rad%n_phy,rad%n_phy,title='Full Hamiltonian Matrix lm-'//itoc(lm))
      END DO
  END IF
!
  Call Diagnalize_Global_Matrix( wa, mat, rad%n_phy )   
END SUBROUTINE Radial_Mat
!**********************************************************************
!***********************************************************************
!deck Theta_Mat
!***begin prologue     Theta_Mat
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Theta_Mat(atom, ang, reg, wa )
  IMPLICIT NONE         
  TYPE(ATOMS)                                      :: atom
  TYPE(THETA)                                      :: ang
  TYPE(REGIONAL), DIMENSION(:)                     :: reg
  TYPE(EVEN_MATRICES)                              :: even
  TYPE(ODD_MATRICES)                               :: odd
  TYPE(FINAL_MATRICES), DIMENSION(:), ALLOCATABLE  :: mat
  TYPE (WORKING_ARRAYS)                            :: wa
  INTEGER                                          :: i
  INTEGER                                          :: lm
  CHARACTER(LEN=3)                                 :: itoc
  CHARACTER(LEN=8)                                 :: coordinate_type
  lower = int_zero
  upper = mtop
  skip = int_two
  DO i = 1,ang%n_reg
     ALLOCATE( reg(i)%mat(lower:upper) )
  END DO
  Call Type_Matrices(reg, even, ang%n_reg)
  Call Normalized_Matrices(reg, even, ang%n_reg)
  Call Final_matrix_Elements( even=even, reg=reg, mat=mat, n_reg=ang%n_reg, type='theta' )
  IF ( mtop > 0 ) THEN
       pre_factor = -one
       Call Type_Matrices(reg, odd, ang%n_reg)
       Call Normalized_Matrices(reg, odd, ang%n_reg)
       lower=int_one
       Call Final_matrix_Elements( odd=odd, reg=reg, mat=mat, n_reg=ang%n_reg, type='theta' )
  END IF
  lower = int_zero
  skip = int_one
  Call Transform_Matrix_to_Standard_Form ( reg, mat, ang%n_reg )
  ALLOCATE( wa%mat(lower:upper) )  
  write(iout,*)
  write(iout,*) '                              Forming the Final Full Matrix'
  write(iout,*)
  Call Form_Global_Matrix(reg, wa, mat, ang%n_reg, ang%n_phy )
  IF (prnt(3) == .true. ) THEN
      DO lm = lower, upper, skip
         Call Print_Matrix(type_real_matrix,wa%mat(lm)%ham,ang%n_phy,ang%n_phy,title='Full Hamiltonian Matrix lm-'//itoc(lm))
      END DO
  END IF
!
  Call Diagonalize_Global_Matrix( wa, mat, ang%n_phy )   
END SUBROUTINE Theta_Mat
!**********************************************************************
!***********************************************************************
!deck Phi_Mat
!***begin prologue     Phi_Mat
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source
!***purpose            Compute sector grid matrix elements.

!***references

!***routines called    iosys, util and mdutil
!***end prologue       

  SUBROUTINE Phi_Mat(atom, ang, reg, wa )
  IMPLICIT NONE         
  TYPE(ATOMS)                                      :: atom
  TYPE(PHI)                                        :: ang
  TYPE(REGIONAL), DIMENSION(:)                     :: reg
  TYPE(EVEN_MATRICES)                              :: even
  TYPE(FINAL_MATRICES), DIMENSION(:), ALLOCATABLE  :: mat
  TYPE (WORKING_ARRAYS)                            :: wa
  INTEGER                                          :: i
  CHARACTER(LEN=3)                                 :: itoc
  CHARACTER(LEN=8)                                 :: coordinate_type
  lower = int_zero
  upper = mtop
  skip = int_two
  DO i = 1,ang%n_reg
     ALLOCATE( reg(i)%mat(lower:upper) )
  END DO
  Call Type_Matrices(reg, even, ang%n_reg)
  Call Normalized_Matrices(reg, even, ang%n_reg)
  Call Final_matrix_Elements( even=even, reg=reg, mat=mat, n_reg=ang%n_reg, type='phi' )
END SUBROUTINE Phi_Mat
!***********************************************************************
!***********************************************************************
!deck Even_Mat
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
!***end prologue       Even_Mat
  SUBROUTINE Even_Mat(reg,even,n_reg)
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
END SUBROUTINE Even_Mat
!***********************************************************************
!***********************************************************************
!deck Odd_Mat
!***begin prologue     Odd_Mat
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

  SUBROUTINE Odd_Mat(reg,odd,n_reg)
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
END SUBROUTINE Odd_Mat
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
  SUBROUTINE Transform_Matrix_to_Standard_Form(reg, mat, n_reg )
  IMPLICIT NONE
  TYPE(REGIONAL), DIMENSION(:)                   :: reg
  TYPE(FINAL_MATRICES), DIMENSION(:)             :: mat
  INTEGER                                        :: n_reg
  INTEGER                                        :: lm
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
  SUBROUTINE Form_Global_Matrix(reg, wa, mat, n_reg, n )
  IMPLICIT NONE
  TYPE (WORKING_ARRAYS)                   :: wa
  TYPE(REGIONAL), DIMENSION(:)            :: reg
  TYPE(FINAL_MATRICES), DIMENSION(:)      :: mat
  INTEGER                                 :: n_reg
  INTEGER                                 :: n
  INTEGER                                 :: lm
  INTEGER                                 :: i
  INTEGER                                 :: first
  INTEGER                                 :: last
!  
  DO lm = lower, upper, skip
     ALLOCATE ( wa%mat(lm)%ham(1:n,1:n) )
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
  SUBROUTINE Diagonalize_Global_Matrix( wa, mat, n )
  IMPLICIT NONE
  TYPE(FINAL_MATRICES), DIMENSION(:)             :: mat
  TYPE(WORKING_ARRAYS)                           :: wa
  INTEGER                                        :: count
  INTEGER                                        :: n
  INTEGER                                        :: i
  INTEGER                                        :: j
  INTEGER                                        :: lm
  INTEGER                                        :: info
  CHARACTER(LEN=3)                               :: itoc
!  
  write(iout,*)
  write(iout,*) '                              Diagonalizing the Full Matrix'
  write(iout,*)
  ALLOCATE(wa%eigenvectors(1:n,1:n),                &
           wa%eigenvalues(1:n),                     &
           wa%work(3*n),                            &
           wa%lower_mat( n*(n+1)/2 ) )
  DO lm = lower, upper, skip
     count = 0
     DO i = 1, n
        DO j = 1, i
           count = count + 1
           wa%lower_mat(count) = wa%mat(lm)%ham(i,j)
        END DO
     END DO
!
!    Used packed form of diagonalizer
!
     Call dspev('v','u',n,wa%lower_mat,wa%eigenvalues,wa%eigenvectors,n,wa%work,info)
     IF (prnt(10) == .true. ) THEN
         Call Print_Matrix(type_real_vector,wa%eigenvalues,title='Eigenvalues lm-'//itoc(lm))
     END IF
     IF (prnt(11) == .true. ) THEN
         Call Print_Matrix(type_real_matrix,wa%eigenvectors,n,n,title='Eigenvectors lm-'//itoc(lm))
     END IF
  END DO
  DEALLOCATE(wa%eigenvectors, wa%eigenvalues, wa%work, wa%lower_mat )
END SUBROUTINE Diagonalize_Global_Matrix
!***********************************************************************
           END MODULE Matrix_Elements
!***********************************************************************
