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
  Subroutine CC_Prop_Input(hamiltonian_source,ham_type,preconditioner,type_storage)
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  USE Channel_Module
  USE CC_Prop_Module
  USE Pack_Global
  USE Iterative_Global
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  CHARACTER (LEN=*)                        :: hamiltonian_source
  CHARACTER (LEN=*)                        :: preconditioner
  CHARACTER (LEN=*)                        :: ham_type
  CHARACTER (LEN=*)                        :: type_storage
  LOGICAL                                  :: dollar
  LOGICAL                                  :: logkey
  LOGICAL                                  :: propagation_test
  REAL*8                                   :: fpkey
  INTEGER                                  :: input, output
  INTEGER                                  :: intkey, i, j
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
  spdim=1
  units='atomic-units'
  hbar=1.d0
  mass=1.d0
  diag=.true.
  typke='dvr'
  coord(1) = 'x'
  IF ( dollar('$cc_setup',card,cpass,inp) ) then
!
       system=chrkey(card,'coordinate_system','cartesian',' ')
       hamiltonian_source=chrkey(card,'hamiltonian_source','internal',' ')
       type_calculation=chrkey(card,'type_calculation','imaginary_time',' ')
       ham_type=chrkey(card,'hamiltonian_type','real',' ')
       type_storage=chrkey(card,'type_storage','full','')
       genpts=logkey(card,'automate_points',.false.,' ')
       proj=logkey(card,'projections',.false.,' ')
       no_cc=logkey(card,'no_coupled_channels',.false.,' ')
       no_pot=logkey(card,'no_one_body_potential',.false.,' ')
       diag_mod=chrkey(card,'diagonal_modification','none',' ')
       non_orth=logkey(card,'non_orthogonal_basis',.false.,' ')
       packed=logkey(card,'pack_matrices',.false.,' ')
       preconditioner=chrkey(card,'preconditioner','none',' ')
       print_parameter=logkey(card,'print=cc=on',.false.,' ')
       print_packed=logkey(card,'print=packed_matrix=on',.false.,' ')
       IF(packed) THEN
          in_core=logkey(card,'in_core',.false.,' ')
          drop_tol = fpkey(card,'drop_tolerance',1.d-08,' ')
          lenbuf = intkey(card,'buffer_length',lenbuf,' ')
          type_packing = chrkey(card,'type_packing','by_column',' ')
          ALLOCATE(mat_var(4))
       END IF
       IF(type_calculation == 'imaginary_time') THEN
          type_calculation='imaginary_time'
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
!
!
!
       con=fpkey(card,'eigenvalue_convergence',1.d-08,' ')
       WRITE(iout,1) type_calculation
!
!       Print options are set here.  Taking the default(none) will give minimal
!       print which is ok most of the time but not when debugging or getting
!       a feel for how the code runs.  Using the all-details option gives lots
!       of information.
!
       pr_main(1)=chrkey(card,'print','none',' ')
       log_main(1)=.false.
       IF(pr_main(1)=='all_details') THEN
          log_main(1)=.true.
       END IF

!
!
!
  ELSE
       write(iout,2)
       stop
  END IF  
!
!
!  Look for the keyword $time in the input file and read all variables
!  that are associated with this keyword.
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
      t_init=fpkey(card,'initial_time',0.D0,' ')
      ntreg=intkey(card,'number_of_time_regions',1,' ')
      deltat=fpkey(card,'time_interval',.01D0,' ')
      nvec=intkey(card,'number_of_initial_wave_packets',1,' ')
      e_method=chrkey(card,'eigenvalue_method','hamiltonian',' ')
      WRITE(iout,3) ntreg, deltat
  ELSE
      write(iout,2)
      stop
  END IF
!
!
  nc=1
  IF(.not.no_cc) THEN
     IF( dollar('$channel',card,cpass,inp) ) then
         nc=intkey(card,'number_of_channels',1,' ')
         Max_L=intkey(card,'L_max',0,' ')
         Max_M=intkey(card,'M_max',0,' ')
         WRITE(iout,4) nc, Max_L, Max_M
         WRITE(iout,5)
         ALLOCATE(channel(nc))
         DO i=1,nc
            ALLOCATE(channel(i)%labels(4))
            CALL intarr(card,'channel_labels_channel_'//itoc(i),channel(i)%labels,4,' ')
            WRITE(iout,6) channel(i)%labels(1:4)
         END DO
     ELSE
         write(iout,2)
         stop
     END IF
  END IF
!
!
1 FORMAT(/,20X,'solve CC_Prop time-dependent schrodinger equation in'    &
         /,20x,'          ',a16)
2 FORMAT(/,5x,'no basis card section')
3 FORMAT(/,5X,'number of time intervals = ',i4,/,5X,                     &
              'time step                = ',5X,F15.8 )
4 FORMAT(/,20X,'Channel Information'                                     &
         /,5x, 'Number of Channels = ',i5,5x,'Maximum L   = ', i3,5x,    &
               'Maximum M = ', i3)
5 FORMAT(/,25x,'Channel Information',//,10x,'   L   ','   S   ',         &
                                            '   M_L  ','   M_S  ')
6 FORMAT(/,9x,i5,2x,i5,3x,i5,3x,i5)
END Subroutine CC_Prop_Input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
