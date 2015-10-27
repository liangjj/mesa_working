!deck CC_Prop_Input.f
!***begin prologue     CC_Prop_Input
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           input, propagation
!***author             schneider, b. i.(nsf)
!***source
!***purpose            input for close-coupled time propagation
!***
!***references
!***routines called    iosys, util and mdutil
!***modules used
!**                    Name                 Location
!                      ----                 -------
!
!                   dvrprop_global          Modules
!                   dvr_shared              Modules
!                   dvr_global              Modules
!
!***end prologue       input
!
!
  Subroutine CC_Prop_Input(hamiltonian_source,ham_type,preconditioner,type_storage,file_directory)
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  USE Pack_Global
  USE Hamiltonian_Module
  USE Lanczos_Global,       only:maximum_number_of_time_subintervals
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  CHARACTER (LEN=*)                        :: hamiltonian_source
  CHARACTER (LEN=*)                        :: preconditioner
  CHARACTER (LEN=*)                        :: ham_type
  CHARACTER (LEN=*)                        :: type_storage
  CHARACTER (LEN=*)                        :: file_directory
  LOGICAL                                  :: dollar
  LOGICAL                                  :: logkey
  LOGICAL                                  :: propagation_test
  REAL*8                                   :: fpkey
  INTEGER                                  :: input, output
  INTEGER                                  :: len
  INTEGER                                  :: intkey
  INTEGER                                  :: i
  COMMON /io/ input, output
!
!  Get the input and output file numbers which appear in the 
!  Module input_output.f90 and put them into input and output which appears
!  in the ONLY common block in the code.  This is needed to pass into library
!  routines and for no other purpose.
!
  input = inp
  output = iout
!
!  Open the input and output files
!
  OPEN(input,file='CC_Prop.inp',status='old')
  OPEN(output,file='CC_Prop.out',status='unknown')
  WRITE(iout,*)
!
!  Look for the keyword $CC_Setup in the input file and read all variables
!  that are associated with this keyword.
!
!***********************************************************************************************
! Main input into setting up the calculation
!***********************************************************************************************
  units='atomic_units'
  hbar=1.d0
  mass=1.d0
  IF ( dollar('$cc_setup',card,cpass,inp) ) then
!
       file_directory=chrkey(card,'file_directory','FILE_DIRECTORY',' ')
       system=chrkey(card,'coordinate_system','cartesian',' ')
       hamiltonian_source=chrkey(card,'hamiltonian_source','internal',' ')
       con=fpkey(card,'eigenvalue_convergence',1.d-08,' ')
!
!***********************************************************************************************
!      Automatically set matrices to real in using imaginary time or complex for
!      real time
!
       type_calculation=chrkey(card,'type_calculation','imaginary_time',' ')
!***********************************************************************************************
!       
       ham_type=chrkey(card,'hamiltonian_type','real',' ')
       IF (type_calculation == 'real_time') THEN
           ham_type = 'complex'
       END IF
       WRITE(iout,1) type_calculation
       IF(type_calculation == 'imaginary_time') THEN
          absorb=.false.
       ELSE
          type_calculation='real_time'
!
!         We might need an absorbing potential
!
          absorb=logkey(card,'add_absorbing_potential',.false.,' ') 
          IF(absorb) THEN
             n_reg_absorb(1)=intkey(card,'number_of_absorbing_'//        &
                                         'regions',0,' ')
          END IF
       END IF
!***********************************************************************************************
       type_storage=chrkey(card,'type_storage','triangle','')
       proj=logkey(card,'projections',.false.,' ')
       no_cc=logkey(card,'no_coupled_channels',.false.,' ')
       no_pot=logkey(card,'no_one_body_potential',.false.,' ')
       diag_mod=chrkey(card,'diagonal_modification','none',' ')
       non_orth=logkey(card,'non_orthogonal_basis',.false.,' ')
       packed=logkey(card,'pack_matrices',.false.,' ')
       preconditioner=chrkey(card,'preconditioner','none',' ')
       IF(.not.non_orth) THEN
          preconditioner = 'none'
       END IF
!***********************************************************************************************
       print_packed=logkey(card,'print_packed_matrix',.false.,' ')
!      If using packed forms, this needs to be set here
!***********************************************************************************************
       IF(packed) THEN
          in_core=logkey(card,'in_core',.false.,' ')
          drop_overlap = fpkey(card,'overlap_drop_tolerance',1.d-08,' ')
          drop_cholesky = fpkey(card,'cholesky_drop_tolerance',1.d-08,' ')
          drop_hamiltonian = fpkey(card,'hamiltonian_drop_tolerance',1.d-08,' ')
          lenbuf = intkey(card,'buffer_length',lenbuf,' ')
          type_packing = chrkey(card,'type_packing','by_column',' ')
          ALLOCATE(mat_var(4))
       END IF
!***********************************************************************************************
!
!       Print options are set here.  Taking the default(none) will give minimal
!       print which is ok most of the time but not when debugging or getting
!       a feel for how the code runs.  Using the all-details option gives lots
!       of information.
!
!***********************************************************************************************
       pr_main(1)=chrkey(card,'print','none',' ')
       log_main(1)=.false.
       IF(pr_main(1)=='all_details') THEN
          log_main(1)=.true.
       END IF
!
  ELSE
       write(iout,2)
       stop
  END IF  
!
!***********************************************************************************************
  write(iout,3) system, hamiltonian_source, type_calculation, ham_type
  IF(packed) THEN
     Write(iout,4) packed, drop_overlap, drop_cholesky, drop_hamiltonian,         &
                   lenbuf, type_packing
  END IF
  IF(preconditioner /= 'none') THEN
     Write(iout,5) preconditioner
  END IF
!***********************************************************************************************
!   Time information for propagation
!***********************************************************************************************
!  Look for the keyword $time in the input file and read all variables
!  that are associated with this keyword.
!***********************************************************************************************
!
  IF( dollar('$time',card,cpass,inp) ) then
!  
!     Read in the number of time regions and delta t and the print options.
!     
      DO  i=3,8
          pr_main(i)='print=main='//pr_main(i)
      END DO
      pr_main(9)=chrkey(card,'print=main=',pr_main(9),' ')
      IF(pr_main(9) == 'all') THEN
         CALL setprn(log_main(3),6)
      ELSE
         CALL setlog(log_main(3),pr_main(3),card,6)
      END IF
!***********************************************************************************************
!     
      t_init=fpkey(card,'initial_time',0.D0,' ')
      ntreg=intkey(card,'number_of_time_regions',1,' ')
      maximum_number_of_time_subintervals = intkey(card,'maximum_number_of_time_subintervals',10,' ')
      deltat=fpkey(card,'time_interval',.01D0,' ')
      nvec=intkey(card,'number_of_initial_wave_packets',1,' ')
      e_method=chrkey(card,'eigenvalue_method','hamiltonian',' ')
      WRITE(iout,6) ntreg, deltat
  ELSE
      write(iout,2)
      stop
  END IF
!***********************************************************************************************
!
!
!
1 FORMAT(/,20X,'Solve CC_Prop Time-Dependent Schrodinger Equation in'            &
         /,20x,'          ',a16)
2 FORMAT(/,5x,'no basis card section')
3 FORMAT(/,5X,'Coordinate System          = ',a16,                               &
           5x,'Hamiltonian Source         = ',a8,                                &
         /,5x,'Type Calculation           = ',a16,                               &
           5x,'Hamiltonian Type           = ',a8)
4 FORMAT(/,5x,'Pack Matrices              = ',l1,14x,                            &   
           5x,'Overlap Drop Tolerance     = ',e15.8,                             &     
         /,5x,'Cholesky Drop Tolerance    = ',e15.8,                             &     
           5x,'Hamiltonian Drop Tolerance = ',e15.8,                             &     
         /,5x,'Buffer Length              = ',i10,5x,                            &
           5x,'Type Packing               = ',a16 )
5 FORMAT(/,5x,'Preconditioner             = ',a32)
6 FORMAT(/,5X,'Number of Time Intervals   = ',i4,                                &
           5X 'Time Step = ',F15.8 )
END Subroutine CC_Prop_Input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
