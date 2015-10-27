!deck SO_driver
!**begin prologue     SO_driver
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           **author             schneider, barry (nsf)
!**source
!**purpose            driver for second order differential equations
!**description        
!**                   
!**references
!**routines called      Name                    Location    
!                       ----                    --------
!                     Input_Prop             Prop_Sources
!                     Space_Prop             Prop_Sources
!                     Regional_Matrices      Prop_Modules:regional_module
!
!**modules used         Name                    Location               
!                       ----                    --------
!                    dvrprop_global          Modules   
!                    Regional_Matrices       Prop_Modules
!
!**end prologue       SO_driver
  PROGRAM SO_driver
  USE dvrprop_global
  USE dvr_global
  USE dvr_shared
  USE regional_module
  USE Coulomb_Functions_Module,    ONLY : type, asymptotic_size, angular_momentum,     &
                                          asymptotic_series, print_sigma_l,            &
                                          print_long_range_coefficients,               &
                                          print_convergence
  IMPLICIT NONE
!
  REAL*8                                   :: dscale = -.5d0
  REAL*8                                   :: fpkey
  INTEGER                                  :: i 
  INTEGER                                  :: j 
  INTEGER                                  :: nen 
  INTEGER                                  :: l_minimum 
  INTEGER                                  :: l_maximum  
  INTEGER                                  :: intkey 
  INTEGER                                  :: iostat 
  INTEGER                                  :: input
  INTEGER                                  :: output
  INTEGER                                  :: Number_of_Data_Sets
  LOGICAL                                  :: Compute_Regional_Matrices
  LOGICAL                                  :: logkey
  LOGICAL                                  :: dollar
  CHARACTER(LEN=80)                        :: chrkey
  CHARACTER(LEN=3)                         :: itoc
  INTEGER                                  :: Number_of_Energies
  REAL*8, DIMENSION(:), ALLOCATABLE        :: energies
  REAL*8, DIMENSION(:), ALLOCATABLE        :: v_hold
  COMMON /io/ input, output
!
  input = inp
  output = iout
!
!  Open the input and output files
!
  OPEN(input,file='so.inp',status='old')
  OPEN(output,file='so.out',status='unknown')
!
! Read in the Input
!
  WRITE(iout,1)
  WRITE(iout,2)
  WRITE(iout,1)
  spdim = 1
  system = 'cartesian'
  coord(1) = 'x'
!
! Calculate the Global FEDVR Matrices
!
  IF ( dollar('$begin_data',card,cpass,inp) ) then
       Number_of_Data_Sets = intkey(card,'number_of_data_sets',1,' ')
       Compute_Regional_Matrices=logkey(card,'compute_regional_matrices',.false.,' ')
  END IF
!
!      Get all the matrices.  We will set up the Hamiltonian ourselves.
!      All we need is the KE matrix at this point.
!
  DO i = 1, Number_of_Data_Sets
     IF ( dollar('$data_set_'//itoc(i),card,cpass,inp) ) then
          charge = fpkey(card,'charge',-1.d0,' ')
          Number_of_Energies = intkey(card,'number_of_energies',1,' ')
          l_minimum = intkey(card,'smallest_angular_momentum',0,' ')
          l_maximum = intkey(card,'largest_angular_momentum',0,' ')
          type = chrkey(card,'equation_type','coulomb',' ')
          keywrd = chrkey(card,'key_word','dvr',' ')
          Number_of_Energies = intkey(card,'number_of_energies',1,' ')
          ALLOCATE(energies(Number_of_Energies))
          Call fparr(card,'energies',energies,Number_of_Energies,' ')
          IF(charge /= zero) THEN
             print_sigma_l = logkey(card,'print_sigma_l',.false.,' ')
             print_long_range_coefficients = logkey(card,                      &
                                           'print_long_range_coefficients',    &
                                            .false.,' ')
             print_convergence = logkey(card,'print_convergence',.false.,' ')
             asymptotic_size = intkey(card,                                    &
                               'number_of_terms_in_asymptotic_series',50,' ')
             ALLOCATE( asymptotic_series(0:l_maximum) )
          END IF
!
!         Compute all the arrays needed which are independent
!         of angular momentum
!
     ALLOCATE(grid(1), num_reg(1), nfun_reg(maxreg,1) )
     CALL dvr_input(nphy(1),nglobal(1),coord(1))
     ALLOCATE(grid(1)%pt(nphy(1)),                             &
              grid(1)%wt(nphy(1)),                             &
              grid(1)%f(nphy(1),nphy(1)),                      &
              grid(1)%df(nphy(1),nphy(1)),                     &   
              grid(1)%ddf(nphy(1),nphy(1)),                    &   
              grid(1)%ke(nphy(1),nphy(1)),                     &   
              grid(1)%p_mom(nphy(1),nphy(1)),                  &   
              grid(1)%h(nphy(1),nphy(1)),                      &   
              grid(1)%v(nphy(1)),grid(1)%srf_prm(2))
     ALLOCATE(grid(1)%eigv_0(nphy(1)),                         &   
              grid(1)%eigvec_0(nphy(1),nphy(1)),               &   
              grid(1)%eigv(nphy(1)),                           &   
              grid(1)%eigvec(nphy(1),nphy(1)),                 &   
              grid(1)%srf_0(nphy(1),2),                        &
              grid(1)%srf(nphy(1),2))    
     ALLOCATE(v_hold(nphy(1)))
     CALL Space_DVR
     row(1) = 2*nphy(1) - 2
     num_reg(1) = nreg
     IF(bcl == 0) then
        npt(1) = npt(1) - 1
     END IF
     IF(bcr == 0 ) then
        npt(nreg) = npt(nreg) - 1
     END IF
     nfun_reg(:,i) = npt
     n_reg_real = num_reg
     key='FEDVR'
     WRITE(iout,3) spdim, key, nphy(1)
     IF (Compute_Regional_Matrices) THEN
         CALL Regional_Matrices
     END IF
!
!          Calculate the potential
!
     Call v_mat(v_hold,grid(1)%pt,i,coord,nphy(1))
!
!               Loop over the angular momentum
!
     DO j=l_minimum,l_maximum
        Write(iout,1)
        angular_momentum = j
        Write(iout,4) grid(1)%pt(nphy(1)), angular_momentum
        grid(1)%v = v_hold          
        call addang(grid(1)%v,grid(1)%pt,dscale,j,coord,       &
                    nphy(1),.false.)
        Call Make_Hamiltonian(grid(1)%ke,grid(1)%v,            &
                              grid(1)%eigv,grid(1)%eigvec,     &
                              grid(1)%srf_prm,                 &
                              grid(1)%srf(:,2),nphy(1))    
        DO nen=1, Number_Of_Energies
           Call Phase_Shift(grid(1)%srf(:,2),grid(1)%eigv,     &
                            grid(1)%pt(nphy(1)),energies(nen), &
                            nphy(1))
        END DO
        Write(iout,1)
     END DO
     DEALLOCATE(grid(1)%pt,                                    &
                grid(1)%wt,                                    &
                grid(1)%f,                                     &
                grid(1)%df,                                    &   
                grid(1)%ddf,                                   &   
                grid(1)%ke,                                    &   
                grid(1)%p_mom,                                 &   
                grid(1)%h,                                     &   
                grid(1)%v,                                     &
                grid(1)%srf_prm)
     DEALLOCATE(grid(1)%eigv_0,                                &   
                grid(1)%eigvec_0,                              &   
                grid(1)%eigv,                                  &   
                grid(1)%eigvec,                                &   
                grid(1)%srf_0,                                 &
                grid(1)%srf)    
     DEALLOCATE(grid, num_reg, nfun_reg)
     IF (Compute_Regional_Matrices ) THEN
         DO j=1,n_reg_real(1)
            DEALLOCATE( mat_reg(j,1)%pt_d, mat_reg(j,1)%ke_mat_d )
         END DO
         DEALLOCATE(mat_reg)
     END IF
          DEALLOCATE(energies)
     END IF
  END DO
!
1 FORMAT('************************************************************************')
2 FORMAT(25X,'Solve Second Order DE')
3 Format(/,1x,'Number of Dimensions   = ',i2,   &
           1x,'Spatial discretization = ',a8,   &
         /,1x,'Number of Points       = ',i10)
4 FORMAT(25X,'Results: R_Matrix Radius = 'e15.8,2x,'Angular Momentum = ',i3)
  stop
END PROGRAM SO_driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
