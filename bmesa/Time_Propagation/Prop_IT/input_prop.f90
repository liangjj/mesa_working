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
!***modules            dvrprop_global_it, dvr_shared, dvr_global
!***end prologue       input
  Subroutine input_prop
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  LOGICAL                                  :: dollar, logkey
  REAL*8                                   :: fpkey
  INTEGER                                  :: intkey, i, j
  CHARACTER (LEN=4), DIMENSION(10)         :: ans 
!
  ans(1) = 'no'
  IF (ans(1) == 'yes') then
      WRITE(5,*) 'propagation method'
      READ(5,*) algorithm
      WRITE(5,*) '     coordinate system, number of space variables,'
      WRITE(5,*) '     and number of time regions'
      READ(5,*) system, spdim, ntreg
      WRITE(5,*) '     coordinate labels for each space variable'
      READ(5,*)  ( coord(i), i=1,spdim )
      WRITE(5,*) '     answer yes or no'
      WRITE(5,*) '     atomic units, automate points, non-linear potentials,' 
      WRITE(5,*) '     no spatial Hamiltonian, diagonalize, imaginary time'
      WRITE(5,*) '     plots and projections'
      READ(5,*) (ans(i),i=1,8)
      IF(ans(1) == 'yes') then
         units='atomic-units'
         hbar=1.d0
         mass=1.d0 
      END IF
      IF(ans(2) == 'yes') then
         genpts=.true.
      ELSE
         genpts=.false.
      END IF
      IF(ans(3) == 'yes') then
         nlse=.true.
      ELSE
         nlse=.false.
      END IF
      IF(ans(4) == 'yes') then
         space=.true.
      ELSE
         space=.false.
      END IF
      IF(ans(5) == 'yes') then
         diag=.true.
      ELSE
         diag=.false.
      END IF
      IF(ans(6) == 'yes') then
         plot=.true.
      ELSE
         plot=.false. 
      END IF
      IF(ans(7) == 'yes') then
         proj=.true.
      ELSE
         proj=.false. 
      END IF
      IF(algorithm /= 'arnoldi') THEN
         WRITE(5,*) '     Sector Print Information'
         READ(5,*) ans(1)
         IF(ans(1) == 'yes') then
            pr_main(1) = 'all-details'
            log_main(1) = .true.
         ELSE
            log_main(1) = .false.
         END IF
         WRITE(5,*) '     Propagation Order'
         READ(5,*) prop_order
         WRITE(5,*) '     Diagonal Modification'
         READ(5,*) diag_mod
      END IF
      WRITE(5,*) '     DVR(dvr or packed) or FD(fd) representation'
      READ(5,*) typke
      plot_step = 1
      IF(plot == .true.) then
         WRITE(5,*) '   Plot Step'
         READ(5,*) plot_step
         WRITE(5,*) plot_step
      END IF   
      WRITE(5,*) 'coordinate system            = ',system
      WRITE(5,*) 'number of space variables    = ', spdim
      WRITE(5,*) 'number of time regions       = ',ntreg
      WRITE(5,*) 'coordinate labels            = ',(coord(i)(1:4),i=1,spdim)
      WRITE(5,*) 'units = ',units(1:16),'automate points = ',genpts
      WRITE(5,*) 'non-linear potentials = ',nlse,'no spatial Hamiltonian = ',space
      WRITE(5,*) 'diagonalize = ',diag 
      WRITE(5,*) 'plots = ',plot,'projections = ',proj 
      IF(algorithm /= 'arnoldi') THEN
         WRITE(5,*) 'sector-print = ',log_main(1)
         WRITE(5,*) 'propagation order = ',prop_order 
         IF (prop_order > 4 ) THEN
             write(iout,1)
             call lnkerr('propagation order error')
         END IF
         WRITE(5,*) 'diagonal modification = ',diag_mod 
      END IF
      WRITE(5,*) 'space representation = ',typke 
  ELSE  
      WRITE(iout,*)
      IF ( dollar('$prop_basis',card,cpass,inp) ) then
!
!         Set spatial dimensionality of problem and program 
!         options
!
           algorithm=chrkey(card,'propagation-method',            &
                                 'split_operator',' ')
           con=fpkey(card,'eigenvalue-convergence',1.d-06,' ')
           spdim=intkey(card,'number-of-space-variables',1,' ')
           WRITE(iout,2) algorithm, spdim
           system=chrkey(card,'coordinate-system','cartesian',' ')
           units=chrkey(card,'units','atomic-units',' ')
           IF(units == 'atomic-units') then
              hbar=1.d0
              mass=1.d0
           END IF
           genpts=logkey(card,'automate-points',.false.,' ')
           nlse=logkey(card,'non-linear-equations',.false.,' ')
           plot=logkey(card,'plot',.false.,' ')
           plot_step=intkey(card,'plot_step',1,' ')
           space=logkey(card,'no-spatial-hamiltonian',.false.,' ')
           diag=logkey(card,'get-eigenpairs',.false.,' ')
           typke=chrkey(card,'kinetic-energy-type','dvr',' ')
           proj=logkey(card,'projections',.false.,' ')
           write(iout,3) system, units, nlse, space, diag, typke, &
                         genpts, proj
           IF(algorithm /= 'arnoldi') THEN
              pr_main(1)=chrkey(card,'print','none',' ')
              log_main(1)=.false.
              IF(pr_main(1)=='all-details') THEN
                 log_main(1)=.true.
              END IF
              diag_mod=chrkey(card,'diagonal-modification','none',' ')
              prop_order=intkey(card,'propagation-order',2,' ')
              IF (prop_order > 4 ) THEN
                  write(iout,1)
                  call lnkerr('propagation order error')
              END IF
              WRITE(iout,4) prop_order, diag_mod
           END IF
           DO i=1, spdim
              coord(i)=chrkey(card,'space-variable-'//itoc(i),'x',' ')
           END DO   
      ELSE
           write(iout,5)
           stop
      END IF  
  END IF
!
  IF( dollar('$time',card,cpass,inp) ) then
!  
!     Read in the number of time regions and delta t
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
      e_method=chrkey(card,'eigenvalue-method','exponential',' ')
      WRITE(iout,6) ntreg, deltat
  END IF

!
!       open BEC file 
!
  call iosys ('read character "bec filename" from rwf',    &
              -1,0,0,filbec)
  call iosys ('open bec as new',0,0,0,filbec)
  call iosys('rewind all on bec read-and-write',0,0,0,' ')
!
1 FORMAT(/,5x,'error in propagtion order')
2 FORMAT(/,20X,'solve time-dependent schrodinger equation in'    &
         /,20x,'          imaginary time',//,20x,                &
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
