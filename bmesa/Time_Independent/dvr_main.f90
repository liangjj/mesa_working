  PROGRAM dvr_main
  USE dvr_shared
  USE dvr_global
  USE dvr_prnt
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  CHARACTER (LEN=80)                       :: chrkey
  LOGICAL                                  :: dollar, logkey, itdiag
  REAL*8                                   :: fpkey, tau
  REAL*4                                   :: secnds
  INTEGER                                  :: intkey, i, j, iostat, nroot
  INTEGER                                  :: n_r
  CHARACTER*4096                           :: ops 
  CHARACTER*4, DIMENSION(10)               :: ans 
  CALL Drum
  CALL IOsys('read character options from rwf',-1,0,0,ops)
  WRITE(iout,1) spdim
  WRITE(iout,2)
  WRITE(5,*) '          Code For Time Independent Schroedinger Equation'
  WRITE(5,*) '          Using DVR or FD Representation'
  WRITE(5,*) '          Do You Wish Terminal Input of Variables'
  READ(5,*) ans(1)
  IF (ans(1) == 'yes') then
      WRITE(5,*) '     coordinate system, number of space dimensions,'
      READ(5,*) system, spdim
      WRITE(5,*) '     coordinate labels for each dimension'
      READ(5,*)  ( coord(i), i=1,spdim )
      WRITE(5,*) '     answer yes or no'
      WRITE(5,*) '     atomic units, automate points, non-linear potentials,' 
      WRITE(5,*) '     iterative diagonalization' 
      READ(5,*) (ans(i),i=1,4)
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
         itdiag = .true.
      ELSE
         itdiag = .false.
      END IF
      WRITE(5,*) '     Sector Print Information'
      READ(5,*) ans(1)
      IF(ans(1) == 'yes') then
         prloc(1) = 'sector-details'
         prn(1) = .true.
      ELSE
         prn(1) = .false.
      END IF
      WRITE(5,*) '     DVR(dvr or packed) or FD(fd) representation'
      READ(5,*) typke
      WRITE(5,*) 'coordinate system = ',system
      WRITE(5,*) 'number of space dimensions = ', spdim
      WRITE(5,*) 'coordinate labels = ',(coord(i)(1:3),i=1,spdim)
      WRITE(5,*) 'units = ',units(1:16),'automate points = ',genpts
      WRITE(5,*) 'non-linear potentials = ',nlse 
      WRITE(5,*) 'sector-print = ',prn(1)
      WRITE(5,*) 'space representation = ',typke 
  ELSE  
      WRITE(iout,*)
      IF ( dollar('$dvr_basis',card,cpass,inp) ) then
!
!         Set spatial dimensionality of problem and program 
!         options
!
           spdim=intkey(card,'number-of-space-dimensions',1,' ')
           genpts=logkey(card,'automate-points',.false.,' ')
           nlse=logkey(card,'non-linear-equations',.false.,' ')
           units=chrkey(card,'units','atomic-units',' ')
           IF(units == 'atomic-units') then
              hbar=1.d0
              mass=1.d0
           END IF
           system=chrkey(card,'coordinate-system','cartesian',' ')
           typke=chrkey(card,'kinetic-energy-type','dvr',' ')
           prloc(1)=chrkey(card,'print','sector-details',' ')
           IF(prloc(1) == 'sector-details') THEN
              prn(1)=.true.
           ELSE
              CALL setlog(prn,prloc,card,1)
           END IF
           DO i=1, spdim
              coord(i)=chrkey(card,'dimension-'//itoc(i),'x',' ')
           END DO   
      ELSE
           write(iout,3)
           stop
      END IF  
  END IF
  diag=.true.
!
!
!       Get all of the one-dimensional matrices needed to construct
!       the spatial part of the hamiltonian and associated quantities.
!
!
! Allocate storage for the element matrices and potential
!   
  tim(1) = secnds(0.0)
  ALLOCATE( grid(spdim), num_reg(spdim), nfun_reg(maxreg,spdim) )
  DO  i=1,spdim
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
         ALLOCATE(grid(i)%eigv_0(nphy(i)),           &    
                  grid(i)%eigvec_0(nphy(i),nphy(i)), &    
                  grid(i)%eigv(nphy(i)),             &    
                  grid(i)%eigvec(nphy(i),nphy(i)),   &    
                  grid(i)%srf_0(nphy(i),2),          &
                  grid(i)%srf(nphy(i),2))    
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
         row(i)=2*nphy(i) - 2
!
!        Deallocate everything except the points, weights, functions,
!        kinetic energy and potential arrays
!
         DEALLOCATE(grid(i)%df,grid(i)%ddf,grid(i)%p_mom,grid(i)%h, &
                    grid(i)%srf_prm)
         DEALLOCATE(grid(i)%eigv,grid(i)%eigvec,grid(i)%srf_0,grid(i)%srf)
         num_reg(i) = nreg
         IF(bcl == 0) then
            npt(1) = npt(1) - 1
         END IF
         IF(bcr == 0 ) then
            npt(nreg) = npt(nreg) - 1
         END IF
         nfun_reg(:,i) = npt
         write(iout,4) i
         write(iout,5) 
         DO j=1,num_reg(i)
            write(iout,6) j, nfun_reg(j,i)
         END DO
      ELSE
!
!           The current code is only set up for a three point FD
!           formula.  
!
         CALL fd_input(nphy(i),nglobal(i),row(i),coord(i))
         ALLOCATE(grid(i)%pt(nphy(i)),               &
                  grid(i)%wt(nphy(i)),               &
                  grid(i)%ke(row(i),nphy(i)),        &
                  grid(i)%h(row(i),nphy(i)),         &
                  grid(i)%v(nphy(i)))
         ALLOCATE(grid(i)%eigv_0(nphy(i)),           &
                  grid(i)%eigvec_0(nphy(i),nphy(i)), &
                  grid(i)%eigv(nphy(i)),             &
                  grid(i)%eigvec(nphy(i),nphy(i)))
         CALL fd_basis(pt_0(i),grid(i)%pt,grid(i)%wt,              &
                       grid(i)%ke,grid(i)%eigv_0,grid(i)%eigvec_0, &
                       grid(i)%h,grid(i)%eigv,grid(i)%eigvec,      &
                       grid(i)%v,nphy(i),nglobal(i),row(i),coord(i))
         DEALLOCATE(grid(i)%h,grid(i)%eigv,grid(i)%eigvec)
!
         IF (row(i) == 2 ) then
             nreg = nphy(i) - 1
             n_r = 2
         ELSE IF(row(i) == 3 ) then
             nreg = nphy(i) - 2
             n_r = 3
         ELSE IF(row(i) == 4 ) then
             nreg = nphy(i) - 3
             n_r = 4
         END IF
         num_reg(i) = nreg
         nfun_reg(:,i) = n_r
         write(iout,7) i
         write(iout,5) 
         DO j=1,num_reg(i)
            write(iout,6) j, nfun_reg(j,i)
         END DO
      END IF
  END DO
  tim(2)=secnds(0.0)
  del(1)=tim(2)-tim(1)
  WRITE(iout,8) del(1)
!
! Get potential parameters.  The actual potential itself
! is calculated in the propagation step.
!
  CALL v_couple
!
! Now do the actual diagonalization
!
  IF(.not. itdiag) then
      nroot=intkey(ops,'number-of-roots',5,' ')
  ELSE
      CALL dvd_dat
  END IF

!  CALL full_diag
!
! Deallocate all arrays that have been allocated for the entire
! calculation.
! 
  DO  i=1,spdim
      DEALLOCATE(grid(i)%pt,grid(i)%wt,grid(i)%ke,grid(i)%v)
      DEALLOCATE(grid(i)%eigv_0,grid(i)%eigvec_0)
  END DO
  call chainx(0)
  stop
1    FORMAT(/,20X,'time-independent basis function code',//,20X,  &
                  'number of spatial dimensions = ',i1)
2    FORMAT(/,15X,'calculation = solve time-independent schrodinger'  &
                  ' equation')
3    FORMAT(/,5x,'no basis card section')
4    FORMAT(/15x,'Regional Information for DVR Basis,' &
                 ' Spatial Dimension = ',i2)
5    FORMAT(/,10x,'Region',5x,'Number of Functions')
6    FORMAT(11x,i2,14x,i3)
7    FORMAT(/15x,'Regional Information for FD Basis,' &
                 ' Spatial Dimension = ',i2)
8    FORMAT('***********************************************' &
            '*************************'                       &
            /,10X,'time to compute the spatial Hamiltonian '  &
            /,10x,'and associated quantities = ',f15.8,/,     &
            '***********************************************' &
            '*************************')
END PROGRAM dvr_main

