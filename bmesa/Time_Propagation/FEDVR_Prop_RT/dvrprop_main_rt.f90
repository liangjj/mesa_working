  PROGRAM dvrprop_main_rt
  USE dvrprop_global_rt
  USE dvr_shared
  USE dvr_global
  USE regional_diag
  USE regional_umat
  USE finite_element_propagator
  USE propagator
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  LOGICAL                                  :: dollar, logkey
  REAL*8                                   :: fpkey
  REAL*4                                   :: secnds
  INTEGER                                  :: intkey, i, j, iostat
  INTEGER                                  :: start, last
  INTEGER                                  :: row_dim
  INTEGER, DIMENSION(2)                    :: words
  CHARACTER (LEN=4096)                     :: ops 
  CHARACTER (LEN=4), DIMENSION(10)         :: ans 
  CHARACTER(LEN=8)                         :: mat_typ
!
  CALL Drum
  CALL IOsys('read character options from rwf',-1,0,0,ops)
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
        diag_mod=chrkey(card,'diagonal-modification','none',' ')
        diag=logkey(card,'get-eigenpairs',.false.,' ')
        typke=chrkey(card,'kinetic-energy-type','dvr',' ')
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
        IF(absorb) THEN
           DO i=1, spdim
              n_reg_absorb(i)=intkey(card,'number-of-absorbing-'//        &
                                       'regions-space-variable-'//     &
                                        itoc(i),0,' ')
           END DO   
        END IF
        IF(nlse) THEN
           CAll vnlse
        END IF
  ELSE
       write(iout,1)
       stop
  END IF  
  IF (prop_order > 4 ) THEN
      write(iout,2)
      call lnkerr('propagation order error')
  END IF
  OPEN (UNIT=iplot(1),FILE='sector_mat',                               &
        ACCESS='sequential',FORM='formatted',                          &
        IOSTAT=IOSTAT,STATUS='unknown')
  IF(IOSTAT /= 0) THEN
     CALL lnkerr('error in file handling')
  END IF
  WRITE(iout,3) spdim
  WRITE(iout,4)
  WRITE(iout,5) prop_order, space, diag_mod, absorb
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
      WRITE(iout,6) ntreg, deltat
  END IF
  time(1)=secnds(0.0)
!
! Allocate storage for the element matrices and potential
!   
  ALLOCATE(grid(spdim), mat_reg_d(maxreg,spdim),                       &
           num_reg(spdim), nfun_reg(maxreg,spdim) )
!
! For scattering problems, we allow for an absorbing potential in the
! last region for each coordinate.
!
  IF(absorb) THEN
     ALLOCATE(mat_reg_z(maxreg,spdim))
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
            call dvr_basis(pt_0(i),                                    &
                           grid(i)%pt,grid(i)%wt,                      &
                           grid(i)%f,grid(i)%df,                       &
                           grid(i)%ddf,                                &
                           grid(i)%ke,                                 &
                           grid(i)%p_mom,                              &
                           grid(i)%eigv_0,                             &
                           grid(i)%eigvec_0,                           &
                           grid(i)%h,                                  &
                           grid(i)%eigv,                               &
                           grid(i)%eigvec,                             &
                           grid(i)%v,                                  &
                           grid(i)%srf_prm,                            &
                           grid(i)%srf_0,                              &
                           grid(i)%srf,                                &
                           coord(i),nphy(i),nglobal(i))
!           Deallocate everything except the points, weights,
!           kinetic energy and potential arrays
!
!
!           Save one body eigenvalues and vectors.
!
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
      n_reg_real(i) = num_reg(i)- n_reg_absorb(i)
      DO j=1,n_reg_real(i)
         write(iout,10) j, nfun_reg(j,i),'real'
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
      IF(absorb) THEN
!
!        There is an absorber.  Allocate the complex arrays.
!
         DO j=n_reg_real(i)+1, num_reg(i)
            write(iout,10) j, nfun_reg(j,i),'complex'
            ALLOCATE(mat_reg_d(j,i)%pt_d                               &
                    (nfun_reg(j,i)),                                   &
                     mat_reg_z(j,i)%v_add_z                            &
                    (nfun_reg(j,i)),                                   &
                     mat_reg_d(j,i)%ke_mat_d                           &
                    (nfun_reg(j,i),nfun_reg(j,i)),                     &
                     mat_reg_z(j,i)%ke_mat_z                           &
                    (nfun_reg(j,i),nfun_reg(j,i)),                     &
                     mat_reg_z(j,i)%eigval_mat_z                       &
                    (nfun_reg(j,i)),                                   &
                     mat_reg_z(j,i)%eigvec_mat_z_r                     &
                    (nfun_reg(j,i),nfun_reg(j,i)),                     &
                     mat_reg_z(j,i)%eigvec_mat_z_l                     & 
                    (nfun_reg(j,i),nfun_reg(j,i)))
            maxmem_reg=max(maxmem_reg,nfun_reg(j,i))
         END DO
         ALLOCATE(scr_z(10*maxmem_reg))
      END IF
      IF(typke == 'dvr'.OR.typke == 'packed') THEN
         mat_typ='full'
         row_dim=nphy(i)
      ELSE
         mat_typ='banded'
         row_dim=row(i)
      END IF
      if(diag_mod /= 'none') THEN
         call modify_diag(grid(i)%ke,grid(i)%v,row_dim,nphy(i),        &
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
      IF(absorb) THEN
         DO j=n_reg_real(i)+1,num_reg(i)
            call add_absorb(mat_reg_d(j,i)%pt_d,                       &
                            mat_reg_d(j,i)%ke_mat_d,                   &
                            mat_reg_z(j,i)%v_add_z,                    &
                            mat_reg_z(j,i)%ke_mat_z,                   &
                            nfun_reg(j,i))
         END DO
      END IF
!
      DO j=1,n_reg_real(i)
         write(iplot(1),*) 'real sector = ', j, 'size = ',npt(j)
         write(iplot(1),*) 'kinetic energy minus diagonals'
         write(iplot(1),*) mat_reg_d(j,i)%ke_mat_d
!
!        Diagonalize
!
         call diag_reg(mat_reg_d(j,i)%ke_mat_d,                        &
                       mat_reg_d(j,i)%eigval_mat_d,                    &
                       mat_reg_d(j,i)%eigvec_mat_d,                    &
                       mat_reg_d(j,i)%eigvec_mat_d,                    &
                       nfun_reg(j,i),j)
      END DO
      DEALLOCATE(scr_d)
      IF(absorb) THEN     
         write(iout,*) n_reg_real(i), num_reg(i)
         DO j=n_reg_real(i)+1,num_reg(i)
            write(iplot(1),*) 'complex sector = ', j, 'size = ',npt(j)
            write(iplot(1),*) 'kinetic energy minus diagonals'
            write(iplot(1),*) mat_reg_z(j,i)%ke_mat_z
            call diag_reg(mat_reg_z(j,i)%ke_mat_z,                     &
                          mat_reg_z(j,i)%eigval_mat_z,                 &
                          mat_reg_z(j,i)%eigvec_mat_z_r,               &
                          mat_reg_z(j,i)%eigvec_mat_z_l,               &
                          nfun_reg(j,i),j)
         END DO  
         DEALLOCATE(scr_z)
      END IF
!        Calculate the time dependent propagators
!
      ALLOCATE(si_d(maxmem_reg,maxmem_reg),                            &
               ci_d(maxmem_reg,maxmem_reg))
      IF(absorb) THEN
         ALLOCATE(si_z(maxmem_reg,maxmem_reg),                         &
                  ci_z(maxmem_reg,maxmem_reg))
      END IF
      p_fac = 1.d0         
      n_prop = 1
      IF(prop_order == 4) THEN
         p_fac = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
         n_prop = 2
      END IF
      DO j=1,num_reg(i)
         ALLOCATE(mat_reg_d(j,i)%cosine_t_mat                          &
                 (nfun_reg(j,i),nfun_reg(j,i),n_prop),                 &
                  mat_reg_d(j,i)%sine_t_mat                            &
                 (nfun_reg(j,i),nfun_reg(j,i),n_prop))
      END DO
      IF(.not.absorb) THEN
        starting_reg = 1
        ending_reg = num_reg(i)
        n_reg = num_reg(i)
        write(iout,20)
        call reg_prop_d(i,deltat)
      ELSE
        starting_reg = 1
        ending_reg = n_reg_real(i)
        n_reg = n_reg_real(i)
        write(iout,20)
        CALL reg_prop_d(i,deltat)
        starting_reg = n_reg_real(i) + 1
        ending_reg = num_reg(i)
        n_reg = ending_reg - starting_reg + 1
        write(iout,30)
        call reg_prop_z(i,.5d0*deltat)
      END IF
      DEALLOCATE(si_d,ci_d)
      IF(absorb) THEN
         DEALLOCATE(si_z,ci_z)
      END IF
      DO j=1,n_reg_real(i)
         DEALLOCATE(mat_reg_d(j,i)%ke_mat_d,                           &
                    mat_reg_d(j,i)%eigvec_mat_d,                       &
                    mat_reg_d(j,i)%eigval_mat_d)
      END DO          
      IF(absorb) THEN
         DO j=n_reg_real(i)+1, num_reg(i)
            DEALLOCATE(mat_reg_z(j,i)%v_add_z,                         &
                       mat_reg_d(j,i)%ke_mat_d,                        &
                       mat_reg_z(j,i)%ke_mat_z,                        &
                       mat_reg_z(j,i)%eigval_mat_z,                    &
                       mat_reg_z(j,i)%eigvec_mat_z_r,                  &
                       mat_reg_z(j,i)%eigvec_mat_z_l)
         END DO
      END IF
  END DO
  time(2)=secnds(0.0)
  delta(1)=time(2)-time(1)
  WRITE(iout,40) delta(1)
!
! Get potential parameters.  The actual potential itself
! is calculated in the propagation step.
!
  CALL v_couple
!
! Now do the actual Propagation
!
! Allocate some needed scratch space.
!
  IF(spdim == 2) THEN
     maxdim=max(nphy(2),nphy(1))
     words(1)=max(nphy(1)*4,maxdim*2)
     ALLOCATE( f_1(nphy(1)), f_2(words(1),1) )
  ELSE IF(spdim == 3) then
     maxdim=max(nphy(3),nphy(2),nphy(1))
     words(1)=max(nphy(2)*nphy(1),maxdim*2,nphy(1)*4)
     words(2)=max(nphy(2)*nphy(1)*4,maxdim*maxdim*2)
     ALLOCATE( f_1(nphy(1)), f_2(words(1),1), f_3(words(2),1,1) )
  END IF
!
! Allocate the main arrays after determining the total dimenisonality.
!
  n3d=1
  DO  i=1,spdim
      n3d=n3d*nphy(i)
  END DO
  ALLOCATE(v_tot(n3d),sin_diag(n3d),cos_diag(n3d),tim_pts(ntreg+1))
  CALL setup
  IF(spdim == 1) THEN
     ALLOCATE(psi_1d(nphy(1),2),v_scr_1d(nphy(1),2))
     CALL so_prop(psi_1d,v_scr_1d)
  ELSE IF(spdim == 2 ) THEN
     ALLOCATE(psi_2d(nphy(2),nphy(1),2),v_scr_2d(nphy(2),nphy(1),2))
     CALL so_prop(psi_2d,v_scr_2d)
  ELSE IF(spdim == 3 ) THEN
     ALLOCATE(psi_3d(nphy(3),nphy(2),nphy(1),2),v_scr_3d(nphy(3),nphy(2),nphy(1),2))
     CALL so_prop(psi_3d,v_scr_3d)
  ELSE
     CALL lnkerr('error in dimension')
  END IF
!
! Deallocate all arrays that have been allocated for the entire
! calculation.
! 
  DEALLOCATE(v_tot,sin_diag,cos_diag,tim_pts)
  IF(spdim == 1) THEN
     DEALLOCATE(psi_1d,v_scr_1d)
  ELSE IF(spdim == 2 ) THEN
     DEALLOCATE(psi_2d,v_scr_2d)
     DEALLOCATE( f_1,f_2 )
  ELSE IF(spdim == 3 ) THEN
     DEALLOCATE(psi_3d,v_scr_3d)
     DEALLOCATE( f_1,f_2,f_3 )
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
        DEALLOCATE(grid(i)%eigv,grid(i)%eigvec,grid(i)%eigv_0,grid(i)%eigvec_0)
     END DO
  END IF 
  DO i=1,spdim
     DO j=1,num_reg(i)
        DEALLOCATE(mat_reg_d(j,i)%cosine_t_mat,                        & 
                      mat_reg_d(j,i)%sine_t_mat)
     END DO
  END DO
  DEALLOCATE( grid,mat_reg_d,num_reg,nfun_reg )
  IF(absorb) THEN
     DEALLOCATE( mat_reg_z )
  END IF
  call chainx(0)
  stop
1    FORMAT(/,5x,'no basis card section')
2    FORMAT(/,5x,'error in propagtion order')
3    FORMAT(/,20X,'time-dependent basis function code',//,20X,         &
                  'number of variables = ',i1)
4    FORMAT(/,15X,'calculation = solve time-dependent schrodinger'     &
                  ' equation')
5    FORMAT(/,25X,'time-dependent data',                               & 
            /,5X,'order of time propagation   = ',i1,                  &
            /,5X,'no spatial hamiltonian      = ',l1,                  &
            /,5x,'diagonal modification       = ',a24,                 &
            /,5x,'absorbing potential added   = ',l1)
6    FORMAT(/,5X,'number of time intervals = ',i4,/,5X,                &
                 'time step                = ',5X,F15.8 )
7    FORMAT(/15x,'Regional Information for DVR Basis,'                 &
                 ' Variable = ',i2)
8    FORMAT(/,10x,'Region',5x,'Number of Functions',5x,'Type')
9    FORMAT(/15x,'Regional Information for FD Basis,'                  &
                 ' Variable = ',i2)
10   FORMAT(11x,i4,14x,i5,7x,a8)
20   FORMAT('Constructing Propagators for Sectors with Real '          &
            'Potentials')
30   FORMAT('Constructing Propagators for Sectors with Complex '       &
            'Potentials')
40   FORMAT('***********************************************'          &
            '*************************'                                &
            /,10X,'time to compute the spatial Hamiltonian '           &
            /,10x,'and associated quantities = ',f15.8,/,              &
            '***********************************************'          &
            '*************************')
END PROGRAM dvrprop_main_rt

