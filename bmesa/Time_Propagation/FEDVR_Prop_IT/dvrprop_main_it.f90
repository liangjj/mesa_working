  PROGRAM dvrprop_main_it
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  USE finite_element_propagator
  USE psi_h_psi
  USE h_on_vector
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  LOGICAL                                  :: dollar, logkey
  REAL*8                                   :: fpkey
  REAL*4                                   :: secnds
  INTEGER                                  :: intkey, i, j, iostat
  INTEGER                                  :: start, last, size
  INTEGER                                  :: row_dim
  INTEGER, DIMENSION(2)                    :: words
  CHARACTER (LEN=4096)                     :: ops 
  CHARACTER (LEN=4), DIMENSION(10)         :: ans 
  CHARACTER(LEN=8)                         :: mat_typ
!
  CALL Drum
  CALL IOsys('read character options from rwf',-1,0,0,ops)
!  WRITE(5,*) '          Code For Imaginary Time-Propagation'
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
         plot=.true.
      ELSE
         plot=.false. 
      END IF
      IF(ans(7) == 'yes') then
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
      WRITE(5,*) '     Propagation Order'
      READ(5,*) prop_order
      WRITE(5,*) '     Diagonal Modification'
      READ(5,*) diag_mod
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
      WRITE(5,*) 'sector-print = ',log_main(1)
      WRITE(5,*) 'propagation order = ',prop_order 
      IF (prop_order > 4 ) THEN
          write(iout,1)
          call lnkerr('propagation order error')
      END IF
      WRITE(5,*) 'diagonal modification = ',diag_mod 
      WRITE(5,*) 'space representation = ',typke 
  ELSE  
      WRITE(iout,*)
      IF ( dollar('$dvrprop_basis',card,cpass,inp) ) then
!
!         Set spatial dimensionality of problem and program 
!         options
!
           spdim=intkey(card,'number-of-space-variables',1,' ')
           WRITE(iout,2) spdim
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
           DO i=1, spdim
              coord(i)=chrkey(card,'space-variable-'//itoc(i),'x',' ')
           END DO   
      ELSE
           write(iout,5)
           stop
      END IF  
  END IF
  IF(nlse) THEN
     CAll vnlse
  END IF
  OPEN (UNIT=iplot(1),FILE='sector_mat',                               &
        ACCESS='sequential',FORM='formatted',                          &
  IOSTAT=IOSTAT,STATUS='unknown')
  IF(IOSTAT /= 0) THEN
     CALL lnkerr('error in file handling')
  END IF
!
!
!       Get all of the one-dimensional matrices needed to construct
!       the spatial part of the hamiltonian and associated quantities.
!
  call iosys ('read character "bec filename" from rwf',                &
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
      nvec=intkey(card,'number-of-initial-wave-packets',1,' ')
      e_method=chrkey(card,'eigenvalue-method','exponential',' ')
      WRITE(iout,6) ntreg, deltat
  END IF
  time(1)=secnds(0.0)
!
! Allocate storage for the element matrices and potential
!   
  ALLOCATE(grid(spdim), mat_reg_d(maxreg,spdim),                       &
           num_reg(spdim), nfun_reg(maxreg,spdim) )
!
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
         IF(diag) then                      
            ALLOCATE(grid(i)%eigv_0(nphy(i)),                          &    
                     grid(i)%eigvec_0(nphy(i),nphy(i)),                &    
                     grid(i)%eigv(nphy(i)),                            &    
                     grid(i)%eigvec(nphy(i),nphy(i)),                  &    
                     grid(i)%srf_0(nphy(i),2),                         &
                     grid(i)%srf(nphy(i),2))    
         END IF
!
         IF(typwt /= 'fourier') THEN
            ALLOCATE(grid(i)%pt(nphy(i)),                              &
                     grid(i)%wt(nphy(i)),                              &
                     grid(i)%f(nphy(i),nphy(i)),                       &
                     grid(i)%df(nphy(i),nphy(i)),                      &   
                     grid(i)%ddf(nphy(i),nphy(i)),                     &   
                     grid(i)%ke(nphy(i),nphy(i)),                      &   
                     grid(i)%p_mom(nphy(i),nphy(i)),                   &   
                     grid(i)%h(nphy(i),nphy(i)),                       &   
                     grid(i)%v(nphy(i)),grid(i)%srf_prm(2))
!
!           Compute the DVR points, weights, functions, first and second
!           derivatives, kinetic energy matrix, eigenvalues and eigenvectors
!           of the kinetic energy matrix, full one-particle Hamiltonian
!           matrix,eigenvalues and eigenvectors of the Hamiltonian matrix, 
!           one-body potential, value of DVR functions at the endpoints, 
!           and value of eigenvectors of the kinetic energy and one-body 
!           Hamiltonian at the endpoints.
!
            call dvr_basis(pt_0(i),                                       &
                           grid(i)%pt,                                    &
                           grid(i)%wt,                                    &
                           grid(i)%f,                                     &
                           grid(i)%df,                                    &
                           grid(i)%ddf,                                   &
                           grid(i)%ke,                                    &
                           grid(i)%p_mom,                                 &
                           grid(i)%eigv_0,                                &
                           grid(i)%eigvec_0,                              &
                           grid(i)%h,                                     &
                           grid(i)%eigv,                                  &
                           grid(i)%eigvec,                                &
                           grid(i)%v,                                     &
                           grid(i)%srf_prm,                               &
                           grid(i)%srf_0,                                 &
                           grid(i)%srf,                                   &
                           coord(i),nphy(i),nglobal(i))
!
!        Deallocate everything except the points, weights, functions,
!        kinetic energy and potential arrays
!
            DEALLOCATE(grid(i)%df,grid(i)%ddf,                         &    
                       grid(i)%p_mom,grid(i)%h,grid(i)%srf_prm)
         ELSE IF(typwt == 'fourier') THEN
            ALLOCATE(grid(i)%pt(nphy(i)),                              &
                     grid(i)%wt(nphy(i)),                              &
                     grid(i)%f(nphy(i),nphy(i)),                       &
                     grid(i)%ke(nphy(i),nphy(i)),                      &   
                     grid(i)%h(nphy(i),nphy(i)),                       &   
                     grid(i)%v(nphy(i)))
            CALL fourier(grid(i)%pt,grid(i)%wt,grid(i)%f,              &
                         edge,nphy(i))
            call fourier_basis(grid(i)%pt,                             &
                               grid(i)%ke,                             &
                               grid(i)%h,                              &
                               grid(i)%v,                              &
                               grid(i)%eigv_0,                         &
                               grid(i)%eigvec_0,                       &
                               grid(i)%eigv,                           &
                               grid(i)%eigvec,                         &
                               coord(i),nphy(i))
!
            DEALLOCATE(grid(i)%h)
         ELSE
            call lnkerr('weight error')
         END IF
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

         IF(diag) then                      
            DEALLOCATE(grid(i)%srf_0,grid(i)%srf)
            call iosys('write real "eigenvalues for variable-'         &
                                    //itoc(i)//'" to bec',nphy(i),     &
                                      grid(i)%eigv,0,' ')
            call iosys('write real "eigenvectors for variable-'        &
                                    //itoc(i)//'" to bec',             &
                                      nphy(i)*nphy(i),                 &
                                      grid(i)%eigvec,0,' ')
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
         ALLOCATE(grid(i)%pt(nphy(i)),                                 &
                  grid(i)%wt(nphy(i)),                                 &
                  grid(i)%ke(row(i),nphy(i)),                          &
                  grid(i)%h(row(i),nphy(i)),                           &
                  grid(i)%v(nphy(i)))
         IF(diag) then
            ALLOCATE(grid(i)%eigv_0(nphy(i)),                          &
                     grid(i)%eigvec_0(nphy(i),nphy(i)),                &
                     grid(i)%eigv(nphy(i)),                            &
                     grid(i)%eigvec(nphy(i),nphy(i)))
         END IF
         CALL fd_basis(pt_0(i),                                        &
                       grid(i)%pt,                                     &
                       grid(i)%wt,                                     &
                       grid(i)%ke,                                     &
                       grid(i)%eigv_0,                                 &
                       grid(i)%eigvec_0,                               &
                       grid(i)%h,                                      &
                       grid(i)%eigv,                                   &
                       grid(i)%eigvec,                                 &
                       grid(i)%v,                                      &
                       nphy(i),nglobal(i),row(i),coord(i))
         DEALLOCATE(grid(i)%h)
         IF(diag) THEN
!
!           Save one body eigenvalues and vectors.
!
          
            call iosys('write real "eigenvalues for variable-'         &
                                    //itoc(i)//'" to bec',nphy(i),     &
                                      grid(i)%eigv,0,' ')
            call iosys('write real "eigenvectors for variable-'        &
                                    //itoc(i)//'" to bec',             &
                                      nphy(i)*nphy(i),                 &
                                      grid(i)%eigvec,0,' ')
!
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
         write(iout,9) i
         write(iout,8) 
      END IF
!
!     Done with the spatial part.
!
      DO j=1,num_reg(i)
         write(iout,10) j, nfun_reg(j,i)
         ALLOCATE(mat_reg_d(j,i)%pt_d                                  &
                  (nfun_reg(j,i)),                                     &
                  mat_reg_d(j,i)%ke_mat_d                              &
                  (nfun_reg(j,i),nfun_reg(j,i)),                       &
                  mat_reg_d(j,i)%eigvec_mat_d                          &
                  (nfun_reg(j,i),nfun_reg(j,i)),                       &
                  mat_reg_d(j,i)%eigval_mat_d                          &
                  (nfun_reg(j,i)))
         maxmem_reg=max(maxmem_reg,nfun_reg(j,i))
      END DO
      ALLOCATE(scr_d(5*maxmem_reg))
      IF(typke == 'dvr'.OR.typke == 'packed') THEN
         mat_typ='full'
         row_dim=nphy(i)
      ELSE
         mat_typ='banded'
         row_dim=row(i)
      END IF
      IF(diag_mod /= 'none') THEN
         call modify_diag(grid(i)%ke,grid(i)%v,row_dim,nphy(i),           &
                          diag_mod,mat_typ)
      END IF
      start = 1
      DO j=1,num_reg(i)
         call pt_reg(grid(i)%pt(start),mat_reg_d(j,i)%pt_d,            &
                     nfun_reg(j,i),nphy(i),j)
         start = start + nfun_reg(j,i) - 1
      END DO
      start = 1
      IF(mat_typ=='full') THEN
         DO j=1,num_reg(i)
            call ke_reg_dvr(grid(i)%ke(start,start),                   &
                            mat_reg_d(j,i)%ke_mat_d,                   &
                            nfun_reg(j,i),nphy(i),j)
            start = start + nfun_reg(j,i) - 1
         END DO
      ELSE IF(mat_typ=='banded') THEN
         DO j=1,num_reg(i)
            call ke_reg_fd(grid(i)%ke(1,start),                        &
                           mat_reg_d(j,i)%ke_mat_d,                    &
                           nfun_reg(j,i),row_dim,j)
            start = start + nfun_reg(j,i) - 1
         END DO
      ELSE
         CALL lnkerr('error')
      END IF  
!
      DO j=1,num_reg(i) 
         write(iplot(1),*) 'sector = ', j, 'size = ',npt(j)
         write(iplot(1),*) 'kinetic energy minus diagonals'
         write(iplot(1),*) mat_reg_d(j,i)%ke_mat_d
!
!        Diagonalize
!
         call diag_reg(mat_reg_d(j,i)%ke_mat_d,                        &
                       mat_reg_d(j,i)%eigval_mat_d,                    &
                       mat_reg_d(j,i)%eigvec_mat_d,                    &
                       nfun_reg(j,i),j)
      END DO
      DEALLOCATE(scr_d)
!
!        Calculate the time dependent propagators
!
      ALLOCATE(exp_tmp_d(maxmem_reg,maxmem_reg))
      p_fac = 1.d0
      n_prop =1
      IF(prop_order == 4) THEN
         p_fac = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
         n_prop = 2
      END IF
      DO j=1,num_reg(i)
         ALLOCATE(mat_reg_d(j,i)%exp_t_mat                            &
                 (nfun_reg(j,i),nfun_reg(j,i),n_prop))
      END DO
      CALL propagator(i)
      DEALLOCATE(exp_tmp_d)
      DO j=1,num_reg(i)
         DEALLOCATE(mat_reg_d(j,i)%ke_mat_d,                           &
                    mat_reg_d(j,i)%eigvec_mat_d,                       &
                    mat_reg_d(j,i)%eigval_mat_d)
      END DO
  END DO
  time(2)=secnds(0.0)
  delta(1)=time(2)-time(1)
  WRITE(iout,20) delta(1)
!
! Get potential parameters.  The actual potential itself
! is calculated in the propagation step.
!
  CALL v_couple
!
! Now do the actual Propagation
!
!
!          Allocate some needed scratch space.
!
  IF(spdim == 2) THEN
     maxdim=max(nphy(2),nphy(1))
     ALLOCATE( f_1(nphy(1)), f_2(maxdim*2,1) )
  ELSE IF(spdim == 3) then
     maxdim=max(nphy(3),nphy(2),nphy(1))
     words(1)=max(nphy(2)*nphy(1),maxdim*2)
     words(2)=maxdim*maxdim*2
     ALLOCATE( f_1(nphy(1)), f_2(words(1),1), f_3(words(2),1,1) )
  END IF
!
! Allocate the main arrays after determining the total dimension
!
  n3d=1
  DO  i=1,spdim
      n3d=n3d*nphy(i)
  END DO
  ALLOCATE(v_tot(n3d),exp_diag(n3d))
  CALL setup
!
!             Allocate the main arrays.
!
!
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     IF (e_method == 'hamiltonian') THEN
!
!      Calculate Energy in propagation using Hamiltonian
!
         ALLOCATE(buf(spdim))
         DO i=1,spdim
            ALLOCATE( buf(i)%d(nphy(i)),              &
                      buf(i)%hbuf(nphy(i)*nphy(i)),   &
                      buf(i)%hibuf(2,nphy(i)*nphy(i)))
            CALL pack_h(i)
         END DO
     END IF
  END IF
  IF(spdim == 1) THEN
     ALLOCATE(psi_1d(nphy(1),nvec),v_scr_1d(nphy(1),nvec))
     CALL so_prop(psi_1d,v_scr_1d)
  ELSE IF(spdim == 2 ) THEN
     ALLOCATE(psi_2d(nphy(2),nphy(1),nvec),v_scr_2d(nphy(2),nphy(1),nvec))
     CALL so_prop(psi_2d,v_scr_2d)
  ELSE IF(spdim == 3 ) THEN
     ALLOCATE(psi_3d(nphy(3),nphy(2),nphy(1),nvec),                    &
              v_scr_3d(nphy(3),nphy(2),nphy(1),nvec))
     CALL so_prop(psi_3d,v_scr_3d)
  ELSE
     CALL lnkerr('error in dimension')
  END IF
!
! Deallocate all arrays that have been allocated for the entire
! calculation.
! 
  DEALLOCATE(exp_diag,tim_pts)
  IF (e_method == 'exponential') THEN
     IF(typke == 'dvr'.OR.typke == 'packed') THEN
         ALLOCATE(buf(spdim))
         DO i=1,spdim
            ALLOCATE( buf(i)%d(nphy(i)), &
                      buf(i)%hbuf(nphy(i)*nphy(i)), &
                      buf(i)%hibuf(2,nphy(i)*nphy(i)))
            CALL pack_h(i)
         END DO
     END IF
     IF(spdim == 1) THEN
        call h_v(psi_1d,v_scr_1d,nvec)
        call check_energy(psi_1d(:,1),v_scr_1d(:,1))
        DEALLOCATE(psi_1d,v_scr_1d)
     ELSE IF(spdim == 2 ) THEN
        call h_v(psi_2d,v_scr_2d,nvec)
        call check_energy(psi_2d(:,:,1),v_scr_2d(:,:,1),f_1)
        DEALLOCATE(psi_2d,v_scr_2d)
        DEALLOCATE( f_1,f_2 )
     ELSE IF(spdim == 3 ) THEN
        call h_v(psi_3d,v_scr_3d,nvec)
        call check_energy(psi_3d(:,:,:,1),v_scr_3d(:,:,:,1),f_1,f_2)
        DEALLOCATE(psi_3d,v_scr_3d)
        DEALLOCATE( f_1,f_2,f_3 )
     END IF
  END IF
  IF(typke /= 'fd') THEN
     DO i=1,spdim
        DEALLOCATE(grid(i)%pt,grid(i)%wt,grid(i)%f,grid(i)%ke,grid(i)%v)
     END DO
  ELSE
     DO i=1,spdim
        DEALLOCATE(grid(i)%pt,grid(i)%wt,grid(i)%ke,grid(i)%v)
     END DO
  END IF
  IF(diag) THEN
     DO i=1,spdim
        DEALLOCATE(grid(i)%eigv,grid(i)%eigvec,grid(i)%eigv_0,    &
                   grid(i)%eigvec_0)
     END DO
  END IF
  DO i=1,spdim
     DO j=1,num_reg(i)
        DEALLOCATE(mat_reg_d(j,i)%exp_t_mat)
     END DO
  END DO
  DEALLOCATE( grid,mat_reg_d,num_reg,nfun_reg )
  IF(typke == 'dvr'.OR.typke == 'packed') THEN
     DO i=1,spdim
        DEALLOCATE( buf(i)%d,                         &
                    buf(i)%hbuf,                      &
                    buf(i)%hibuf)
     END DO
     DEALLOCATE(buf)
  END IF
  call chainx(0)
  stop
1    FORMAT(/,5x,'error in propagtion order')
2    FORMAT(/,20X,'solve time-dependent schrodinger equation in'    &
                  '          imaginary time',//,20X,                                  'number of variables = ',i1)
3    FORMAT(/,25X,'time-dependent data',                            &
            /,5X,'coordinate system          = ',a32,               &
            /,5X,'units                      = ',a32,               &
            /,5X,'non-linear potential       = ',l1,                &
            /,5X,'no spatial hamiltonian     = ',l1,                &
            /,5X,'calculate eigenvalues      = ',l1,                &
            /,5X,'kinetic energy type        = ',a16,               &
            /,5X,'automatic point generation = ',l1,                &
            /,5X,'calculate projections      = ',l1)
4    FORMAT(/,5X,'propagation order          = ',i1,                &
            /,5x,'diagonal modification      = ',a24)
5    FORMAT(/,5x,'no basis card section')
6    FORMAT(/,5X,'number of time intervals = ',i4,/,5X,                &
                 'time step                = ',5X,F15.8 )
7    FORMAT(/15x,'Regional Information for DVR Basis,'                 &
                 ' Variable = ',i2)
8    FORMAT(/,10x,'Region',5x,'Number of Functions')
9    FORMAT(/15x,'Regional Information for FD Basis,'                  &
                 ' Variable = ',i2)
10   FORMAT(11x,i4,14x,i5)
20   FORMAT('***********************************************'          &
            '*************************'                                &
            /,10X,'time to compute the spatial Hamiltonian '           &
            /,10x,'and associated quantities = ',f15.8,/,              &
            '***********************************************'          &
            '*************************')

END PROGRAM dvrprop_main_it

