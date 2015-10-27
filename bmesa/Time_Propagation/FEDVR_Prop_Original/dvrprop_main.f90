  PROGRAM dvrprop_main
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  USE dvr_propagators
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  LOGICAL                                  :: dollar, logkey
  REAL*8                                   :: fpkey, tau
  REAL*4                                   :: secnds
  INTEGER                                  :: intkey, i, j, iostat
  INTEGER                                  :: start, last
  CHARACTER (LEN=4096)                     :: ops 
  CHARACTER (LEN=4), DIMENSION(10)         :: ans 
  CHARACTER(LEN=80), DIMENSION(9)          :: pr_loc 
  CHARACTER(LEN=8)                         :: mat_typ
!
  CALL Drum
  CALL IOsys('read character options from rwf',-1,0,0,ops)
!  WRITE(5,*) '          Code For Time-Propagation'
!  WRITE(5,*) '          Using DVR or FD Representation'
!  WRITE(5,*) '          Do You Wish Terminal Input of Variables'
!  READ(5,*) ans(1)
  ans(1) = 'no'
  IF (ans(1) == 'yes') then
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
         imtime=.true.
      ELSE
         imtime=.false.  
      END IF
      IF(ans(7) == 'yes') then
         plot=.true.
      ELSE
         plot=.false. 
      END IF
      IF(ans(8) == 'yes') then
         proj=.true.
      ELSE
         proj=.false. 
      END IF
      WRITE(5,*) '     Sector Print Information'
      READ(5,*) ans(1)
      IF(ans(1) == 'yes') then
         pr_main(1) = 'all-details'
         log_main(1) = .true.
      ELSE
         log_main(1) = .false.
      END IF
      WRITE(5,*) '     DVR(dvr or packed) or FD(fd) representation'
      READ(5,*) typke
      WRITE(5,*) '     Propagation Order'
      READ(5,*) prop_order
      WRITE(5,*) '     Keep Diagonals'
      READ(5,*) keep_diag
      plot_step = 1
      IF(plot == .true.) then
         WRITE(5,*) '   Plot Step'
         READ(5,*) plot_step
         WRITE(5,*) plot_step
      END IF   
      WRITE(5,*) '     Absorbing Potential'
      READ(5,*) absorb
      WRITE(5,*) 'coordinate system            = ',system
      WRITE(5,*) 'number of space variables    = ', spdim
      WRITE(5,*) 'number of time regions       = ',ntreg
      WRITE(5,*) 'coordinate labels            = ',(coord(i)(1:4),i=1,spdim)
      WRITE(5,*) 'units = ',units(1:16),'automate points = ',genpts
      WRITE(5,*) 'non-linear potentials = ',nlse,'no spatial Hamiltonian = ',space
      WRITE(5,*) 'diagonalize = ',diag,'imaginary time = ',imtime 
      WRITE(5,*) 'plots = ',plot,'projections = ',proj 
      WRITE(5,*) 'sector-print = ',log_main(1)
      WRITE(5,*) 'space representation = ',typke 
      WRITE(5,*) 'propagation order = ',typke 
      WRITE(5,*) 'keep full diagonal = ',keep_diag 
      WRITE(5,*) 'absorbing potential added = ',absorb 
  ELSE  
      WRITE(iout,*)
      IF ( dollar('$dvrprop_basis',card,cpass,inp) ) then
!
!         Set spatial dimensionality of problem and program 
!         options
!
           spdim=intkey(card,'number-of-space-variables',1,' ')
           prop_order=intkey(card,'propagation-order',2,' ')
           genpts=logkey(card,'automate-points',.false.,' ')
           nlse=logkey(card,'non-linear-equations',.false.,' ')
           units=chrkey(card,'units','atomic-units',' ')
           IF(units == 'atomic-units') then
              hbar=1.d0
              mass=1.d0
           END IF
           system=chrkey(card,'coordinate-system','cartesian',' ')
           plot=logkey(card,'plot',.false.,' ')
           plot_step=intkey(card,'plot_step',1,' ')
           space=logkey(card,'no-spatial-hamiltonian',.false.,' ')
           keep_diag=logkey(card,'keep-diagonals',.false.,' ')
           diag=logkey(card,'get-eigenpairs',.false.,' ')
           typke=chrkey(card,'kinetic-energy-type','dvr',' ')
           imtime=logkey(card,'imaginary-time',.false.,' ')
           pr_main(1)=chrkey(card,'print','none',' ')
           absorb=logkey(card,'add-absorbing-potential',.false.,' ')
           log_main(1)=.false.
           IF(pr_main(1)=='all-details') THEN
              log_main(1)=.true.
           END IF
           proj=logkey(card,'projections',.false.,' ')
           DO i=1, spdim
              coord(i)=chrkey(card,'space-variable-'//itoc(i),'x',' ')
           END DO   
      ELSE
           write(iout,5)
           stop
      END IF  
  END IF
  OPEN (UNIT=iplot(1),FILE='sector_mat', &
        ACCESS='sequential',FORM='formatted', &
  IOSTAT=IOSTAT,STATUS='unknown')
  IF(IOSTAT /= 0) THEN
     CALL lnkerr('error in file handling')
  END IF
  WRITE(iout,1) spdim
  WRITE(iout,2)
  WRITE(iout,3) prop_order, space, keep_diag, absorb
!
!
!       Get all of the one-dimensional matrices needed to construct
!       the spatial part of the hamiltonian and associated quantities.
!
  call iosys ('read character "bec filename" from rwf', &
              -1,0,0,filbec)
  call iosys ('open bec as new',0,0,0,filbec)
  call iosys('rewind all on bec read-and-write',0,0,0,' ')
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
      WRITE(iout,6) ntreg, deltat
  END IF
  time(1)=secnds(0.0)
!
! Allocate storage for the element matrices and potential
!   
  ALLOCATE(grid(spdim), mat_reg_d(maxreg,spdim), num_reg(spdim), &
           nfun_reg(maxreg,spdim) )
!
! For scattering problems, we allow for an absorbing potential in the
! last region for each coordinate.
!
  IF(absorb) THEN
     ALLOCATE(mat_reg_z(spdim))
  END IF
  DO  i=1,spdim
!
!-----------------------------------------------------------------------
!
!                   Begin the Code for the DVR Case
!
!-----------------------------------------------------------------------
      IF(typke == 'dvr'.OR.typke == 'packed') THEN
         CALL dvr_input(nphy(i),nglobal(i),coord(i))
!
!         Allocate the needed DVR arrays.
!
!
         ALLOCATE(grid(i)%pt(nphy(i)),            &
                  grid(i)%wt(nphy(i)),            &
                  grid(i)%f(nphy(i),nphy(i)),     &
                  grid(i)%df(nphy(i),nphy(i)),    &    
                  grid(i)%ddf(nphy(i),nphy(i)),   &    
                  grid(i)%ke(nphy(i),nphy(i)),    &   
                  grid(i)%p_mom(nphy(i),nphy(i)), &   
                  grid(i)%h(nphy(i),nphy(i)),     &     
                  grid(i)%v(nphy(i)),grid(i)%srf_prm(2))
         IF(diag) then                      
            ALLOCATE(grid(i)%eigv_0(nphy(i)),           &    
                     grid(i)%eigvec_0(nphy(i),nphy(i)), &    
                     grid(i)%eigv(nphy(i)),             &    
                     grid(i)%eigvec(nphy(i),nphy(i)),   &    
                     grid(i)%srf_0(nphy(i),2),          &
                     grid(i)%srf(nphy(i),2))    
         END IF
!
!           Compute the DVR points, weights, functions, first and second
!           derivatives, kinetic energy matrix, eigenvalues and eigenvectors
!           of the kinetic energy matrix, full one-particle Hamiltonian
!           matrix,eigenvalues and eigenvectors of the Hamiltonian matrix, 
!           one-body potential, value of DVR functions at the endpoints, 
!           and value of eigenvectors of the kinetic energy and one-body 
!           Hamiltonian at the endpoints.
!
         call dvr_basis(pt_0(i),grid(i)%pt,grid(i)%wt,grid(i)%f,grid(i)%df,  &
                        grid(i)%ddf,grid(i)%ke,grid(i)%p_mom,grid(i)%eigv_0, &
                        grid(i)%eigvec_0,grid(i)%h,grid(i)%eigv,             &
                        grid(i)%eigvec,grid(i)%v,grid(i)%srf_prm,            &
                        grid(i)%srf_0,grid(i)%srf,coord(i),nphy(i),nglobal(i))
         write(iplot(1),*) nphy(i)
         write(iplot(1),*) 'total points'
         write(iplot(1),*) grid(i)%pt
         write(iplot(1),*) 'total weights'
         write(iplot(1),*) grid(i)%wt
         write(iplot(1),*) 'total potential'
         write(iplot(1),*) grid(i)%v
         row(i) = 2*nphy(i) - 2
!
!        Save one body eigenvalues and vectors.
!
!
!        Deallocate everything except the points, weights, functions,
!        kinetic energy and potential arrays
!
         DEALLOCATE(grid(i)%df,grid(i)%ddf, &    
                    grid(i)%p_mom,grid(i)%h,grid(i)%srf_prm)
         IF(diag) then                      
            call iosys('write real "eigenvalues for variable-'     &
                                    //itoc(i)//'" to bec',nphy(i), &
                                      grid(i)%eigv,0,' ')
            call iosys('write real "eigenvectors for variable-'    &
                                    //itoc(i)//'" to bec',         &
                                      nphy(i)*nphy(i),             &
                                      grid(i)%eigvec,0,' ')
            DEALLOCATE(grid(i)%eigv,grid(i)%eigvec,grid(i)%srf_0,grid(i)%srf)
         END IF

         num_reg(i) = nreg
         IF(bcl == 0) then
            npt(1) = npt(1) - 1
         END IF
         IF(bcr == 0 ) then
            npt(nreg) = npt(nreg) - 1
         END IF
         nfun_reg(:,i) = npt
         write(iout,7) i
         write(iout,8) 
      ELSE
!-------------------------------------------------------------------------------
!
!                   Begin the Code for the FD Case
!
!------------------------------------------------------------------------------
!           The current code is only set up for a three point FD
!           formula.  
!
!
         CALL fd_input(nphy(i),nglobal(i),row(i),coord(i))
         ALLOCATE(grid(i)%pt(nphy(i)),        &
                  grid(i)%wt(nphy(i)),        &
                  grid(i)%ke(row(i),nphy(i)), &
                  grid(i)%h(row(i),nphy(i)),  &
                  grid(i)%v(nphy(i)))
         IF(diag) then
            ALLOCATE(grid(i)%eigv_0(nphy(i)),           &
                     grid(i)%eigvec_0(nphy(i),nphy(i)), &
                     grid(i)%eigv(nphy(i)),             &
                     grid(i)%eigvec(nphy(i),nphy(i)))
         END IF
         CALL fd_basis(pt_0(i),grid(i)%pt,grid(i)%wt,              &
                       grid(i)%ke,grid(i)%eigv_0,grid(i)%eigvec_0, &
                       grid(i)%h,grid(i)%eigv,grid(i)%eigvec,      &
                       grid(i)%v,nphy(i),nglobal(i),row(i),coord(i))
         DEALLOCATE(grid(i)%h)
         IF(diag) THEN
!
!           Save one body eigenvalues and vectors.
!
          
            call iosys('write real "eigenvalues for variable-'     &
                                    //itoc(i)//'" to bec',nphy(i), &
                                      grid(i)%eigv,0,' ')
            call iosys('write real "eigenvectors for variable-'    &
                                    //itoc(i)//'" to bec',         &
                                      nphy(i)*nphy(i),             &
                                      grid(i)%eigvec,0,' ')
!
            DEALLOCATE(grid(i)%eigv,grid(i)%eigvec)
         END IF
         IF (row(i) == 2 ) then
             nreg = nphy(i) - 1
             maxmem_reg = 2
         ELSE IF(row(i) == 3 ) then
             nreg = nphy(i) - 2
             maxmem_reg = 3
         ELSE IF(row(i) == 4 ) then
             nreg = nphy(i) - 3
             maxmem_reg = 4
         END IF
         num_reg(i) = nreg
         nfun_reg(:,i) = row(i)
         write(iout,10) i
         write(iout,8) 
      END IF
!
!     Done with the spatial part.
!
      DO j=1,num_reg(i)
         write(iout,9) j, nfun_reg(j,i)
         ALLOCATE(mat_reg_d(j,i)%ke_mat_d(nfun_reg(j,i),nfun_reg(j,i)),     &
                  mat_reg_d(j,i)%eigvec_mat_d(nfun_reg(j,i),nfun_reg(j,i)), &
                  mat_reg_d(j,i)%eigval_mat_d(nfun_reg(j,i)))
         maxmem_reg=max(maxmem_reg,nfun_reg(j,i))
      END DO
      IF(absorb) THEN
!
!        There is an absorber.  Allocate the complex arrays.
!
         ALLOCATE(mat_reg_z(i)%ke_mat_z                                      &
                 ( nfun_reg(num_reg(i),i),nfun_reg(num_reg(i),i) ),          &
                  mat_reg_z(i)%eigvec_mat_z_r                                &
                 ( nfun_reg(num_reg(i),i),nfun_reg(num_reg(i),i) ),          &
                  mat_reg_z(i)%eigvec_mat_z_l                                &
                 ( nfun_reg(num_reg(i),i),nfun_reg(num_reg(i),i) ),          &
                  mat_reg_z(i)%eigval_mat_z(nfun_reg(num_reg(i),i)) )
         ALLOCATE(scr_d(5*maxmem_reg),scr_z(10*maxmem_reg))
      ELSE
         ALLOCATE(scr_d(5*maxmem_reg))
      END IF
      IF(typke == 'dvr'.OR.typke == 'packed') THEN
         mat_typ='full'
      ELSE
         mat_typ='banded'
      END IF
      call modify_diag(grid(i)%ke,grid(i)%v,nphy(i),keep_diag,mat_typ)
      call ke_reg(grid(i)%ke,nphy(i),i,mat_typ)
      IF(absorb) THEN
         call add_absorb(mat_reg_d(num_reg(i),i)%ke_mat_d,                   &
                         mat_reg_z(i)%ke_mat_z,                              &
                         nfun_reg(num_reg(i),i))
      END IF
!
      DO j=1,num_reg(i) - 1
         write(iplot(1),*) 'sector = ', j, 'size = ',npt(j)
         write(iplot(1),*) 'kinetic energy minus diagonals'
         write(iplot(1),*) mat_reg_d(j,i)%ke_mat_d
!
!        Diagonalize
!
         call diag_reg_d(mat_reg_d(j,i)%ke_mat_d,                 &
                         mat_reg_d(j,i)%eigval_mat_d,             &
                         mat_reg_d(j,i)%eigvec_mat_d,             &
                         mat_reg_d(j,i)%eigvec_mat_d,             &
                         nfun_reg(j,i),j)
      END DO
      j=num_reg(i) 
      IF(.not.absorb) THEN
         write(iplot(1),*) 'sector = ', j, 'size = ',npt(j)
         write(iplot(1),*) 'kinetic energy minus diagonals'
         write(iplot(1),*) mat_reg_d(j,i)%ke_mat_d
         call diag_reg_d(mat_reg_d(j,i)%ke_mat_d,                 &
                         mat_reg_d(j,i)%eigval_mat_d,             &
                         mat_reg_d(j,i)%eigvec_mat_d,             &
                         mat_reg_d(j,i)%eigvec_mat_d,             &
                         nfun_reg(j,i),j)
         DEALLOCATE(scr_d)
      ELSE
         CALL diag_reg_z(mat_reg_z(i)%ke_mat_z,                   &
                         mat_reg_z(i)%eigval_mat_z,               &
                         mat_reg_z(i)%eigvec_mat_z_r,             &
                         mat_reg_z(i)%eigvec_mat_z_l,             &
                         nfun_reg(j,i),j)
         DEALLOCATE(scr_d,scr_z)
      END IF
!
!        Calculate the time dependent propagators
!
      IF(prop_order ==2) then
         CALL propagator_2(i)
      ELSE IF(prop_order ==4) then
!         CALL propagator_4(i)
      END IF
      DO j=1,num_reg(i)
         DEALLOCATE(mat_reg_d(j,i)%ke_mat_d,                    &
                    mat_reg_d(j,i)%eigvec_mat_d,                &
                    mat_reg_d(j,i)%eigval_mat_d)
      END DO          
      IF(absorb) THEN
         DEALLOCATE(mat_reg_z(i)%ke_mat_z,                      &
                    mat_reg_z(i)%eigval_mat_z,                  &
                    mat_reg_z(i)%eigvec_mat_z_l,                &
                    mat_reg_z(i)%eigvec_mat_z_r)
      END IF
  END DO
  time(2)=secnds(0.0)
  delta(1)=time(2)-time(1)
  WRITE(iout,4) delta(1)
!
! Get potential parameters.  The actual potential itself
! is calculated in the propagation step.
!
  CALL v_couple
!
! Now do the actual Propagation
!
  CALL so_prop
!
! Deallocate all arrays that have been allocated for the entire
! calculation.
! 
  DO i=1,spdim
     DEALLOCATE(grid(i)%pt,grid(i)%wt,grid(i)%ke,grid(i)%v)
     IF(diag) THEN
        DEALLOCATE(grid(i)%eigv_0,grid(i)%eigvec_0)
     END IF
     DO j=1,nreg
        DEALLOCATE(mat_reg_d(j,i)%cosine_t_mat,mat_reg_d(j,i)%sine_t_mat)
     END DO
  END DO
  DEALLOCATE( grid,mat_reg_d,num_reg,nfun_reg )
  IF(absorb) THEN
     DEALLOCATE( mat_reg_z )
  END IF
  call chainx(0)
  stop
1    FORMAT(/,20X,'time-dependent basis function code',//,20X,  &
                  'number of variables = ',i1)
2    FORMAT(/,15X,'calculation = solve time-dependent schrodinger'  &
                  ' equation')
3    FORMAT(/,25X,'time-dependent data',              & 
            /,5X,'order of time propagation   = ',i1, &
            /,5X,'no spatial hamiltonian      = ',l1, &
            /,5x,'keeping full diagonal       = ',l1, &
            /,5x,'absorbing potential added   = ',l1)
4    FORMAT('***********************************************' &
            '*************************'                       &
            /,10X,'time to compute the spatial Hamiltonian '  &
            /,10x,'and associated quantities = ',f15.8,/,     &
            '***********************************************' &
            '*************************')
5    FORMAT(/,5x,'no basis card section')
6    FORMAT(/,5X,'number of time intervals = ',i4,/,5X,  &
                 'time step                = ',5X,F15.8 )
7    FORMAT(/15x,'Regional Information for DVR Basis,' &
                 ' Variable = ',i2)
8    FORMAT(/,10x,'Region',5x,'Number of Functions')
9    FORMAT(11x,i4,14x,i5)
10   FORMAT(/15x,'Regional Information for FD Basis,' &
                 ' Variable = ',i2)
END PROGRAM dvrprop_main

