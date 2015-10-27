!deck input.f
!***begin prologue     input
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           input, propagation
!***author             schneider, b. i.(nsf)
!***source
!***purpose            input for time propagation
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
  Subroutine input_prop
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  LOGICAL                                  :: dollar, logkey
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
!  Open the input and output files when not using the mesa script
!
!  OPEN(input,file='prop.inp',status='old')
!  OPEN(output,file='prop.out',status='unknown')
  WRITE(iout,*)
!
!  Look for the keyword $prop_basis in the input file and read all variables
!  that are associated with this keyword.
!
  IF ( dollar('$prop_basis',card,cpass,inp) ) then
!
!         Set spatial dimensionality of problem and program 
!         options
!
        type_calculation=chrkey(card,'type-calculation','imaginary-time',' ')
        IF(type_calculation == 'imaginary-time') THEN
           type_calculation='imaginary time'
           absorb=.false.
        ELSE
           type_calculation='real time'
!
!          We might need an absorbing potential
!
           absorb=logkey(card,'add-absorbing-potential',.false.,' ') 
        END IF
!
!       Choose an algorithm; split operator or arnoldi
!
        algorithm=chrkey(card,'propagation-method',            &
                              'split_operator',' ')
!
!       Is this a 1D, 2D or 3D calculation.
!
        spdim=intkey(card,'number-of-space-variables',1,' ')
        WRITE(iout,2) type_calculation, algorithm, spdim
        system=chrkey(card,'coordinate-system','cartesian',' ')
!
!       The Module grid_global.f90 defines a whole lot of constants that are
!       used in the calculation.  Are we in atomic units?  If not then the
!       definitions of the constants give hbar in joule-sec and mass of the
!       electron in kilograms.
!
        units=chrkey(card,'units','atomic-units',' ')
        IF(units == 'atomic-units') then
           hbar=1.d0
           mass=1.d0
        END IF
!
!       To automate the way points are generated, set genpts to .true.
!       How this is done is defined elsewhere.
!
        genpts=logkey(card,'automate-points',.false.,' ')
!
!       The parameters below need to be set or the default values
!       will be used.  In most cases the default can be used as there
!       is typically no non-linear potential, the basis is DVR (fd is the
!       other possibility and the convergence criterion is tight enough)
!       
        nlse=logkey(card,'non-linear-equations',.false.,' ')
        space=logkey(card,'no-spatial-hamiltonian',.false.,' ')
        diag=logkey(card,'get-eigenpairs',.false.,' ')
        typke=chrkey(card,'kinetic-energy-type','dvr',' ')
        proj=logkey(card,'projections',.false.,' ')
        con=fpkey(card,'eigenvalue-convergence',1.d-08,' ')
        write(iout,3) system, units, nlse, space, diag, typke, &
                      genpts, proj
!
!       Print options are set here.  Taking the default(none) will give minimal
!       print which is ok most of the time but not when debugging or getting
!       a feel for how the code runs.  Using the all-details option gives lots
!       of information.
!
        IF(algorithm /= 'arnoldi') THEN
           pr_main(1)=chrkey(card,'print','none',' ')
           log_main(1)=.false.
           IF(pr_main(1)=='all-details') THEN
              log_main(1)=.true.
           END IF
!
!          One can modify the diagonal part of the one-body operators to
!          include the diagonal part of the one-body potential.  This is not
!          always important but when you are splitting the operator could result
!          in more accuracy to second order.
!
           diag_mod=chrkey(card,'diagonal-modification','none',' ')
           prop_order=intkey(card,'propagation-order',2,' ')
           IF (prop_order > 4 ) THEN
               write(iout,1)
               call lnkerr('propagation order error')
           END IF
           WRITE(iout,4) prop_order, diag_mod
        END IF
!
!       Now we begin the real work.
!       For each dimension compute the minimal information to proceed.
!
        DO i=1, spdim
!
!          This is just a label to get to a keyword for data entry.
!
           coord(i)=chrkey(card,'space-variable-'//itoc(i),'x',' ')
        END DO   
        IF(absorb) THEN
           DO i=1, spdim
              n_reg_absorb(i)=intkey(card,'number-of-absorbing-'//        &
                                          'regions-space-variable-'//     &
                                           itoc(i),0,' ')
           END DO
        END IF
  ELSE
        write(iout,5)
        stop
  END IF  
!
!  Look for the keyword $time in the input file and read all variables
!  that are associated with this keyword.
!
  IF( dollar('$time',card,cpass,inp) ) then
!  
!     Read in the number of time regions and delta t and the print options.
!     The e_method variable can be used to compute the energy in imaginary
!     time propagations either using the exponential or Hamiltonian form.
!     The latter appears to be more accurate and is the default.
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
      t_init=fpkey(card,'initial-time',0.D0,' ')
      ntreg=intkey(card,'number-of-time-regions',1,' ')
      deltat=fpkey(card,'time-interval',.01D0,' ')
      nvec=intkey(card,'number-of-initial-wave-packets',1,' ')
      e_method=chrkey(card,'eigenvalue-method','hamiltonian',' ')
      WRITE(iout,6) ntreg, deltat
  END IF
!
!
1 FORMAT(/,5x,'error in propagtion order')
2 FORMAT(/,20X,'solve time-dependent schrodinger equation in'    &
         /,20x,'          ',a16,//,20x,                          &
               'algorithm = ',a24,/,20x,                         &
               'number of variables = ',i1)
3 FORMAT(/,25X,'time-dependent data',                            &
         /,5X,'coordinate system          = ',a32,               &
         /,5X,'units                      = ',a32,               &
         /,5X,'non-linear potential       = ',l1,                &
         /,5X,'no spatial hamiltonian     = ',l1,                &
         /,5X,'calculate eigenvalues      = ',l1,                &
         /,5X,'kinetic energy type        = ',a16,               &
         /,5X,'automatic point generation = ',l1,                &
         /,5X,'calculate projections      = ',l1)
4 FORMAT(/,5X,'propagation order          = ',i1,                &
         /,5x,'diagonal modification      = ',a24)
5 FORMAT(/,5x,'no basis card section')
6 FORMAT(/,5X,'number of time intervals = ',i4,/,5X,             &
              'time step                = ',5X,F15.8 )
END Subroutine input_prop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
