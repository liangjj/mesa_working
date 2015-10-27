!***********************************************************************
!**begin prologue     Data_Subroutines_Module
!**date written       080612   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time propagation
!**
!**author             schneider, b. i.(nsf)
!**source             
!**purpose            This module packages all of the routines needed
!***                  to read in information on the matrices and to reformat them.
!***                  It does not perform the reformatting.  That is left for another module.
!***                  
!***                  
!***description       Here are the routines designed to read in the basic data needed to perform
!***                  a time propagation.  A few notes.  The laser-atom interation is described in
!***                  the dipole approximation as 
!***                      V = Sum_i E(t) . r_i
!***                  where the electric field does not depend on the spatial coordinates of the atom.
!***                  This can be expressed in terms of the spherical tensors T_z, T_+ and T_- or
!***                  the cartesian z=T_z, x= -1/sqrt(2) [ T_+ - T_- ], y = i/sqrt(2) [ T_= + T_- ]
!***                  Using the Wigner-Eckert Theorem, the matrix elements between these irreducible
!***                  tensors is,
!***                                                    { L 1 L' }
!***                  < LM | T_q | L'M' > = (-1)**(L-M) {        } < L || T || L'>
!***                                                    { -M q M'}
!***                  where < L || T || L'> is the reduced dipole matrix element and the symbol is curly
!***                  brackets is a 3J symbol.  For linear polarization along z, q = 0 and if the initial
!***                  state has M=0, so does the final state.  This is assumed in the code but it could be changed
!***                  easilt if needed.
!**references
!**modules needed     Input_Output_Module, Input_Data_Module, prop_global_module( only the matrix dimension, n3d )
!**end prologue       Data_Subroutines_Module
!***********************************************************************
!***********************************************************************
!
                           MODULE Data_Subroutines_Module
!
                           USE Atomic_Matrices
                           USE Global_Time_Propagation_Module
                           USE Input_Output
                           USE prop_global,     ONLY: n3d
!
                           IMPLICIT NONE
!
!***********************************************************************
!***********************************************************************
                             CONTAINS
!***********************************************************************
!***********************************************************************
!deck Directives
!***begin prologue     Directives
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Directives
!***
!***references
!***routines called
!***end prologue       Directives
!
  SUBROUTINE Directives(ham_type,file_directory)
  IMPLICIT NONE
  CHARACTER(LEN=*)                      :: ham_type
  CHARACTER(LEN=*)                      :: file_directory
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  LOGICAL                               :: dollar
  LOGICAL                               :: logkey
  INTEGER                               :: intkey
  INTEGER                               :: i
  INTEGER                               :: j
  REAL*8                                :: fpkey
  INTEGER                               :: len
  INTEGER                               :: lenth
  INTEGER                               :: IOSTAT
  IF ( dollar('$directives',card,cpass,inp) ) then
       matrix_directive=chrkey(card,'matrix_directive','by-states',' ')
       ham_file=chrkey(card,'hamiltonian_file_name','from_file',' ')
       lenbuf=intkey(card,'buffer_length',1000000,' ')
       diagonalize_only=logkey(card,'diagonalize_only',.false.,' ')
       full_cholesky_to_disk=logkey(card,'full_cholesky_to_disk',.false.,' ')
       overlap_to_disk=logkey(card,'overlap_to_disk',.false.,' ')
       rows_to_print = intkey(card,'rows_to_print',20,' ')
       eigenvectors_to_print = intkey(card,'eigenvectors_to_print',5,' ')
       print_cc=logkey(card,'print_cc',.false.,' ')
       print_packed_matrices=logkey(card,'print_packed_matrices',.false.,' ')
       print_input_matrices=logkey(card,'print_input_matrices',.false.,' ')
       print_internal_matrices=logkey(card,'print_internal_matrices',.false.,' ')
       print_buffers=logkey(card,'print_buffers',.false.,' ')
       reformatting_control=chrkey(card,'reformatting_control','off',' ')
       reformat_input_matrix_only=logkey(card,'reformat_input_matrix_only',.false.,' ')
       eigenvectors_to_print = intkey(card,'eigenvectors_to_print',5,' ')
       to_standard_eigenvalue_problem = logkey(card,                                       & 
                                          'to_standard_eigenvalue_problem',.false.,' ')
       write(iout,1) matrix_directive, ham_file, packed_matrices, reformatting_control,    &
                     reformat_input_matrix_only, diagonalize_only, full_cholesky_to_disk,  &
                     overlap_to_disk, print_cc, print_packed_matrices,                     &
                     print_input_matrices, print_internal_matrices, print_buffers
  END IF
  IF ( matrix_directive == 'by_channels' ) THEN
       Call read_channel_input(ham_type,file_directory)
  ELSE IF( matrix_directive == 'by_states' ) THEN
       Call read_state_input
  END IF
1 Format (/,5x,'Input Matrices                            = ',a16,                         &
          /,5x,'Hamiltonian File                          = ',a80,                         &
          /,5x,'Packed Matrices                           = ',l1,                          &
          /,5x,'Reformatting Control                      = ',a3,                          &
          /,5x,'Reformat Input Matrix Only                = ',l1,                          &
          /,5x,'Diagonalize Only                          = ',l1,                          &
          /,5x,'Write Full Cholesky Decomposition to Disk = ',l1,                          &
          /,5x,'Write Overlap to Disk                     = ',l1,                          &
          /,5x,'Print Channel Information                 = ',l1,                          &
          /,5x,'Print Packed Matrices                     = ',l1,                          &
          /,5x,'Print input Matrices                      = ',l1,                          &
          /,5x,'Print Internal Matrices                   = ',l1,                          &
          /,5x,'Print Buffers                             = ',l1)
END SUBROUTINE Directives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!***********************************************************************
!deck read_channel_input
!***begin prologue     read_channel_input
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            read input
!***
!***references
!***routines called
!***end prologue       read_input
!
  SUBROUTINE read_channel_input(ham_type,file_directory)
  IMPLICIT NONE
  CHARACTER (LEN=*)                     :: ham_type
  CHARACTER (LEN=*)                     :: file_directory
  CHARACTER(LEN=80)                     :: chrkey
  CHARACTER (LEN=8)                     :: itoc
  LOGICAL                               :: dollar
  LOGICAL                               :: logkey
  INTEGER                               :: intkey
  INTEGER                               :: i
  INTEGER                               :: j
  REAL*8                                :: fpkey
  INTEGER                               :: len
  INTEGER                               :: len_1
  INTEGER                               :: lenth
  INTEGER                               :: IOSTAT
!
!
  len=lenth(ham_file)
  len_1=lenth(file_directory)
  WRITE(iout,*) '    Hamiltonian File Name = '//file_directory(1:len_1)//     &
                                             '/'//ham_file(1:len)
  WRITE(iout,*) '    Non Orthogonal Basis = ',non_orth
  OPEN(UNIT=50,FILE=file_directory(1:len_1)//'/'//ham_file(1:len),            &
               ACCESS='sequential', FORM='unformatted',IOSTAT=IOSTAT,STATUS='old')
  IF(input_matrices(1:6) == 'packed') THEN
!
!         This type of read is for a B-spline matrix where only
!         the non-zero matrix elements have been stored.  There is no channel
!         information on the file.
!
     READ(50) n3d
  ELSE IF(input_matrices == 'oleg_channel') THEN
!
!         This type of read is for a B_spline matrix having explicit
!         information about the channels.  Very easy
!         to get sub-matrices channel by channel.
!
      Call Read_Channel_B_Spline_Data
  END IF
  CLOSE(50)
  WRITE(iout,1) n3d
  tri_size = n3d * ( n3d + 1 ) / 2
  IF(in_core) THEN
     lenbuf = n3d*(n3d-1)/2
  ELSE
     lenbuf = min(lenbuf,n3d*(n3d-1)/2)
  END IF
  write(iout,2) tri_size, lenbuf
  rows_to_print = min(rows_to_print,n3d)      
  eigenvectors_to_print = min(eigenvectors_to_print,n3d)
  WRITE(iout,*) '    Hamiltonian Type    = ',ham_type(1:len)
  IF (packed_matrices) THEN
      WRITE(iout,3)
  END IF
!
1  FORMAT(/,5x,'Input Matrix Size = ',i8)
2  FORMAT(/,5x,'Triangle_Size = ',i8,2x,'Buffer Size = ',i8)
3  Format(/,20x,'Using Packed Forms for all Matrices')
END SUBROUTINE read_channel_input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!***********************************************************************
!deck read_state_input
!***begin prologue     read_state_input
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            read input from B-spline program and reformat matrices
!***
!***description        the input matrices are made up of channel-channel submatrices and
!***                   channel-bound and bound-bound submatrices.  symbolically it looks like
!***                                   H_CC 
!***                                   H_BC H_BB
!***                   only the triangle is stored where the matrices are symmetric such as the Hamiltonian
!***                   and overlap.  the program removes unwanted splines associated with the channel matrices
!***                   which, if left in, would lead to physically incorrect results.  splines at the origin and
!***                   at large distances are excised from the input and the matrices are stored in IOsys format.
!***references
!***routines called
!***end prologue       read_input
!
  SUBROUTINE read_state_input
  IMPLICIT NONE
  CHARACTER(LEN=80)                     :: chrkey
  LOGICAL                               :: dollar
  LOGICAL                               :: logkey
  INTEGER                               :: intkey
!
!
  IF ( dollar('$state_parameters',card,cpass,inp) ) THEN
       L_Max=intkey(card,'L_Max',1,' ')
       number_of_splines=intkey(card,'number_of_splines',500,' ')
       spline_order=intkey(card,'spline_order',8,' ')
       remove_last_spline=logkey(card,'remove_last_spline',.false.,' ')
       remove_next_to_last_spline=logkey(card,'remove_next_to_last_spline',.false.,' ')
       diagonalize_matrix=logkey(card,'diagonalize_matrix',.false.,' ')
       dipole_matrices=logkey(card,'dipole_matrices',.false.,' ')
       write_pointers=logkey(card,'write_pointers',.false.,' ')
       write_channel_labels=logkey(card,'write_channel_labels',.false.,' ')
       Write(iout,1) L_Max, number_of_splines, spline_order, remove_last_spline,   &
                     remove_next_to_last_spline
 
 END IF
1 Format(/,5x,'Reformatting B Spline Matrices:  Maximum L Value            = ',i3, &
         /,5x,'                                 Number of Splines          = ',i4, &
         /,5x,'                                 Spline Order               = ',i2, &
         /,5x,'                                 Remove_Last_Spline         = ',l1, &
         /5x, '                                 Remove_Next_To_Last_Spline = ',l1)
END SUBROUTINE read_state_input
!***********************************************************************
!***********************************************************************
!Deck Read_Channel_B_Spline_Data
!***begin prologue     Read_Channel_B_Spline_Data
!***date written       061201   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Read in B spline data in Oleg's format
!***
!***references
!***routines called
!***end prologue       Read_Channel_B_Spline_Data
!
  SUBROUTINE Read_Channel_B_Spline_Data
  IMPLICIT NONE
  INTEGER                               :: nch_i
  INTEGER                               :: nch_j
  INTEGER                               :: k_cp
  INTEGER                               :: ns
  INTEGER                               :: njunk
  INTEGER, DIMENSION(:), ALLOCATABLE    :: iprm
!
  Read(50) nch_j, nch_i, k_cp
  Write(iout,*) nch_i, nch_j, k_cp
  ALLOCATE(spline_array_i(0:nch_i), spline_array_j(0:nch_j))
  Read(50) spline_array_j, spline_array_i, ns, n3d, n3d
  write(iout,*) spline_array_j, spline_array_i, ns, n3d
  ALLOCATE(iprm(n3d))
  READ(50) iprm
!  write(iout,*) iprm
  DEALLOCATE(spline_array_i, spline_array_j, iprm)
!***********************************************************************
  END SUBROUTINE Read_Channel_B_Spline_Data
!***********************************************************************
!***********************************************************************
!deck Atomic_Prop_Input.f
!***begin prologue     Atomic_Prop_Input
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
  Subroutine Atomic_Prop_Input(ham_type,preconditioner,file_directory)
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  USE Pack_Global
  USE Lanczos_Global,       only:maximum_number_of_time_subintervals
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  CHARACTER (LEN=*)                        :: preconditioner
  CHARACTER (LEN=*)                        :: ham_type
  CHARACTER (LEN=*)                        :: file_directory
  LOGICAL                                  :: dollar
  LOGICAL                                  :: logkey
  REAL*8                                   :: fpkey
  INTEGER                                  :: input, output
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
  OPEN(input,file='Atomic_Prop.inp',status='old')
  OPEN(output,file='Atomic_Prop.out',status='unknown')
  WRITE(iout,*)
!
!  Look for the keyword $CC_Setup in the input file and read all variables
!  that are associated with this keyword.
!
!***********************************************************************************************
!                  Main input into setting up the calculation
!***********************************************************************************************
!
!                Set atomic units as default
!
  units='atomic_units'
  hbar=1.d0
  mass=1.d0
  IF ( dollar('$atomic_prop_setup',card,cpass,inp) ) then
!
       file_directory=chrkey(card,'file_directory','FILE_DIRECTORY',' ')
       cholesky_file_name = chrkey(card,'cholesky_matrix_file_name',      &
                                   'Cholesky_Deomposition',' ')
       state_file_name = chrkey(card,'state_matrix_file_name',            &
                                'State_Matrices',' ')
       con=fpkey(card,'eigenvalue_convergence',1.d-08,' ')
       use_atomic_symmetry=chrkey(card,'use_atomic_symmetry','off',' ')
!
!***********************************************************************************************
!      Automatically set matrices to real in using imaginary time or complex for
!      real time
!
       type_calculation=chrkey(card,'type_calculation','imaginary_time',' ')
!***********************************************************************************************
!       
       ham_type=chrkey(card,'hamiltonian_type','real',' ')
       WRITE(iout,1) type_calculation
!***********************************************************************************************
!
!            More variables that need to be set.  The default is usually  &
!            ok for most of them.
!
       non_orth=logkey(card,'non_orthogonal_basis',.false.,' ')
       input_matrices=chrkey(card,'input_matrices',                       &
                             'packed_in_drake_format',' ')
       output_matrices=chrkey(card,'output_matrices',                     &
                              'packed_in_triangular_iosys_format',' ')
       IF (input_matrices(1:6) == 'packed') THEN
           output_matrices='packed_in'//output_matrices(10:80)
           packed_matrices = .true.
       END IF
       preconditioner=chrkey(card,'preconditioner','none',' ')
       IF(non_orth) THEN
          preconditioner = 'cholesky'
       END IF
!***********************************************************************************************
!      If using packed forms, this needs to be set here
!***********************************************************************************************       
       in_core=logkey(card,'in_core',.false.,' ')
       drop_overlap = fpkey(card,'overlap_drop_tolerance',1.d-08,' ')
       drop_cholesky = fpkey(card,'cholesky_drop_tolerance',1.d-08,' ')
       drop_hamiltonian = fpkey(card,'hamiltonian_drop_tolerance',1.d-08,' ')
       drop_dipole = fpkey(card,'dipole_drop_tolerance',1.d-08,' ')
       drop_vectors = fpkey(card,'vector_drop_tolerance',1.d-08,' ')
       packed_file_name = chrkey(card,'packed_matrix_file_name',          &
                                 'Packed_Matrices',' ')
       ALLOCATE(mat_var(4))
!***********************************************************************************************
!
!       Some print options are set here.  Taking the default(none) will   &
!       give minimal print which is ok most of the time but not when      &
!       debugging or getting a feel for how the code runs.  Using the     &
!        all-details option gives lotsof information.
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
  write(iout,3) type_calculation, ham_type
  Write(iout,4) drop_overlap, drop_cholesky, drop_hamiltonian, drop_dipole,            &
                   input_matrices, output_matrices
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
      time_dependent_potential = logkey(card,'time_dependent_potential',  &
                                        .false.,' ')
      photon_energy = fpkey(card,'photon_energy_eV',au_in_eV,' ')
      photon_energy = photon_energy / au_in_eV      
      time_one_optical_cycle = two_pi/photon_energy
      number_optical_cycles=intkey(card,'number_of_optical_cycles',10,' ')
      steps_per_optical_cycle=intkey(                                     &
                                     card,'number_of_steps_per_'//        &
                                          'optical_cycle',10,' ')
      ntreg = number_optical_cycles * steps_per_optical_cycle
      deltat = time_one_optical_cycle/ntreg 
      t_init=fpkey(card,'initial_time',0.D0,' ')
      maximum_number_of_time_subintervals = intkey(                       &
                                                   card,                  &
                                                   'maximum_number_of_'// &
                                                   'time_subintervals',10,' ')
      nvec=intkey(card,'number_of_initial_wave_packets',1,' ')
      e_method=chrkey(card,'eigenvalue_method','hamiltonian',' ')
      peak_type=chrkey(card,'peak_type','intensity',' ')
      IF (peak_type == 'field' ) THEN
          peak_electric_field=fpkey(card,'peak_electric_field',           &
                                    peak_electric_field,' ')
          peak_intensity = peak_electric_field * peak_electric_field      &
                                              * electric_field_to_intensity
      ELSE
          peak_intensity = fpkey(card,'peak_intensity',                   &
                                       peak_intensity,' ')
          peak_electric_field = sqrt ( peak_intensity /                   &
                                       electric_field_to_intensity )
      END IF
      pulse_duration = fpkey(card,'pulse_duration',pulse_duration,' ')
      ramp_time = fpkey(card,'ramp_time',ramp_time,' ')
      WRITE(iout,6) t_init, ntreg, deltat, photon_energy,                 &
                    time_one_optical_cycle, peak_electric_field,          &
                    peak_intensity, pulse_duration, ramp_time
  ELSE
      write(iout,2)
      stop
  END IF
!***********************************************************************************************
!
!
!
1 FORMAT(/,20X,'Solve Atomic_Prop Time-Dependent Schrodinger Equation in'        &
         /,20x,'          ',a16)
2 FORMAT(/,5x,'no basis card section')
3 FORMAT(/,5x,'Type Calculation           = ',a16,                               &
           5x,'Hamiltonian Type           = ',a8)
4 FORMAT(/,5x,'Overlap Drop Tolerance     = ',e15.8,                             &     
           5x,'Cholesky Drop Tolerance    = ',e15.8,                             &     
         /,5x,'Hamiltonian Drop Tolerance = ',e15.8,                             &     
           5x,'Dipole Drop Tolerance      = ',e15.8,                             &     
         /,5x,'Input Matrix Format        = ',a64,                               &
         /,5x,'Output Matrix Format       = ',a64 )
5 FORMAT(/,5x,'Preconditioner             = ',a32)
6 FORMAT(/,5X,'Initial Time  = ',F15.8,                                          &
           5X,'Number of Time Intervals  = ',i4,                                 &
         /,5X,'Time Step(au) = ',F15.8,                                          &
           5X,'Photon Energy(au)        = ',F15.8,                               &
         /,5X,'Time for One Optical Cycle(au) = ',F15.8,                         &
         /,5X,'Peak Electric Field(au)        = ',F15.8,                         &
         /,5X,'Peak Laser Intensity(W/cm**2)  = ',E15.8,                         &
         /,5x,'Pulse Duration(optical cycles) = ',F15.8,                         &
         /,5x,'Ramp Time(optical cycles)      = ',F15.8)                           
END Subroutine Atomic_Prop_Input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!***********************************************************************
  END  MODULE Data_Subroutines_Module
!***********************************************************************
!***********************************************************************
