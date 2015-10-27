!deck driver.f
!***begin prologue     driver
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           lobatto functions
!***author             schneider, b. i.(nsf)
!***source
!***purpose            1. calculate piecewise lobatto dvr functions and
!***                      their one-body matrices,
!***                   2. calculate the regional, one-body kinetic energy
!***                      matrices.
!***                   3. diagonalize one-body kinetic energy matrices

!***references

!***routines called    iosys, util and mdutil
!***end prologue       driver
  PROGRAM driver
!
! It is important to look at the modules.  They contain many of the variables
! that are used over and over again, constants and most importantly, the allocated
! arrays and pointer variables.  Type variables which are of the pointer variety
! are used extensively in the code.  They are linked and allocated here but 
! defined in the modules.
!
  USE dvr_global
  USE dvr_shared
  USE dvrprop_global_rt
  IMPLICIT NONE
  INTEGER                                 :: intkey, i, j, k, err
  INTEGER                                 :: inp, iout 
  CHARACTER (LEN=2)                       :: itoc
  CHARACTER (LEN=80)                      :: chrkey
  CHARACTER (LEN=8)                       :: mat_typ
  LOGICAL                                 :: dollar, logkey
  INTEGER                                 :: n_alloc, row_dim, start
!
!  This common is put in because old F77 library routines used it
!
  COMMON /io/ inp, iout
  inp=input
  iout=output
!
  open(input,file='dvr.inp',status='old')
  open(output,file='dvr.out',status='unknown')
!
! Keyword input is employed using features developed in the Mesa code
! to manipulate character variables.
! The format is free and the user does not have to worry about order.
!
  IF( dollar ('$dvr_input',card,cpass,input ) ) then
      spdim=intkey(card,'number-of-space-variables',1,' ')
      typke=chrkey(card,'kinetic-energy-type','dvr',' ')
      diag_mod=chrkey(card,'diagonal-modification','none',' ')
      diag=logkey(card,'get-eigenpairs',.false.,' ')
      typke=chrkey(card,'kinetic-energy-type','dvr',' ')
      pr_main(1)=chrkey(card,'print','none',' ')
      log_main(1)=.false.
      IF(pr_main(1)=='all-details') THEN
         log_main(1)=.true.
      END IF
!
      ALLOCATE(grid(spdim),mat_reg_d(maxreg,spdim),                     &
               num_reg(spdim), nfun_reg(maxreg,spdim) )
!
!     Loop over space dimensions.
!
      DO i=1,spdim
         coord(i)=chrkey(card,'space-variable-'//itoc(i),'x',' ')      
         WRITE(output,1) coord(i)
         IF(typke == 'dvr'.OR.typke == 'packed') THEN
!
!           DVR function input routine
!
            call dvr_input(nphy(i),nglobal(i),coord(i))
            WRITE(output,2) nphy(i), nglobal(i)
!
!           Allocation of eigenvalues and eigenvectors
!
            IF(diag) then
                ALLOCATE(grid(i)%eigv_0(nphy(i)),grid(i)%eigv(nphy(i)), &
                         grid(i)%eigvec_0(nphy(i),nphy(i)),             &
                         grid(i)%eigvec(nphy(i),nphy(i)),               & 
                         grid(i)%srf_prm(2),                            &
                         grid(i)%srf_0(nphy(i),2),                      &
                         grid(i)%srf(nphy(i),2),                        &
                         stat=err)
                IF (err /= 0 ) then
                    call lnkerr('allocation error')
                END IF
            END IF
!
!           Allocation of points, weights, kinetic energy etc. matrix elements.
!
            IF(typwt /= 'fourier') THEN
                ALLOCATE(grid(i)%pt(nphy(i)),                           &
                         grid(i)%wt(nphy(i)),                           &
                         grid(i)%f(nphy(i),nphy(i)),                    &
                         grid(i)%df(nphy(i),nphy(i)),                   &
                         grid(i)%ddf(nphy(i),nphy(i)),                  &
                         grid(i)%ke(nphy(i),nphy(i)),                   &
                         grid(i)%p_mom(nphy(i),nphy(i)),                &
                         grid(i)%h(nphy(i),nphy(i)),                    &
                         grid(i)%v(nphy(i)),                            &
                         stat=err)
                IF (err /= 0 ) then
                    call lnkerr('allocation error')
                END IF
!              
!               Calculation of all of the DVR quantities
!
                call dvr_basis(pt_0(i),grid(i)%pt,grid(i)%wt,           &
                               grid(i)%f,grid(i)%df,grid(i)%ddf,        &
                               grid(i)%ke,grid(i)%p_mom,                &
                               grid(i)%eigv_0,grid(i)%eigvec_0,         &
                               grid(i)%h,grid(i)%eigv,                  &
                               grid(i)%eigvec,grid(i)%v,                &
                               grid(i)%srf_prm,grid(i)%srf_0,           &
                               grid(i)%srf,coord(i),nphy(i),nglobal(i))
!
!               Some deallocations for unneeded things
!
                DEALLOCATE(grid(i)%df,                                  &
                           grid(i)%ddf,                                 &
                           grid(i)%p_mom,                               &
                           grid(i)%h)

               IF(diag) then
                   DEALLOCATE(grid(i)%srf_prm,                          & 
                              grid(i)%srf_0,                            &
                              grid(i)%srf)
               END IF
            ELSE
!
!               ALlocations for Fourier basis
!
                ALLOCATE(grid(i)%pt(nphy(i)),                           &
                         grid(i)%wt(nphy(i)),                           &
                         grid(i)%f(nphy(i),nphy(i)),                    &
                         grid(i)%ke(nphy(i),nphy(i)),                   &
                         grid(i)%h(nphy(i),nphy(i)),                    &
                         grid(i)%v(nphy(i)),                            &
                         stat=err)
                IF (err /= 0 ) then
                    call lnkerr('allocation error')
                END IF
                call fourier(grid(i)%pt,grid(i)%wt,                     &
                             grid(i)%f,edge,nphy(i))
                call fourier_basis(grid(i)%pt,                          &
                                   grid(i)%ke,                          &
                                   grid(i)%h,                           &
                                   grid(i)%v,                           &
                                   grid(i)%eigv_0,                      &
                                   grid(i)%eigvec_0,                    &
                                   grid(i)%eigv,                        &
                                   grid(i)%eigvec,                      &
                                   coord(i),nphy(i))
!
!               Deallocations for Fourier basis
!
                DEALLOCATE(grid(i)%pt,                                  &
                           grid(i)%wt,                                  &
                           grid(i)%f,                                   &
                           grid(i)%ke,                                  &
                           grid(i)%h,                                   &
                           grid(i)%v)
            END IF
            IF(diag) then
               DEALLOCATE(grid(i)%eigv_0,                               &
                          grid(i)%eigv,                                 &
                          grid(i)%eigvec_0,                             &
                          grid(i)%eigvec)
            END IF
            row(i) = 2*nphy(i) - 2
            num_reg(i) = nreg
            IF(bcl == 0) then
               npt(1) = npt(1) - 1
            END IF
            IF(bcr == 0 ) then
               npt(nreg) = npt(nreg) - 1
            END IF
            nfun_reg(:,i) = npt
            WRITE(output,4) i
            WRITE(output,5)
         ELSE
!
!              Can also do finite difference type calculations but
!              these have not been as extensively tested.
!
!              Similar philosophy to DVR
!
               CALL fd_input(nphy(i),nglobal(i),row(i),coord(i))
               ALLOCATE(grid(i)%pt(nphy(i)),                            &
                        grid(i)%wt(nphy(i)),                            &
                        grid(i)%ke(row(i),nphy(i)),                     &
                        grid(i)%h(row(i),nphy(i)),                      &
                        grid(i)%v(nphy(i)))
               IF(diag) then
                  ALLOCATE(grid(i)%eigv_0(nphy(i)),                     &
                           grid(i)%eigvec_0(nphy(i),nphy(i)),           &
                           grid(i)%eigv(nphy(i)),                       &
                           grid(i)%eigvec(nphy(i),nphy(i)))
               END IF
               CALL fd_basis(pt_0(i),grid(i)%pt,grid(i)%wt,             &
                             grid(i)%ke,grid(i)%eigv_0,                 &
                             grid(i)%eigvec_0,                          &
                             grid(i)%h,grid(i)%eigv,                    &
                             grid(i)%eigvec, &
                             grid(i)%v,nphy(i),nglobal(i),row(i),coord(i))
               DEALLOCATE(grid(i)%h)
               IF(diag) then
                  DEALLOCATE(grid(i)%eigv_0,                            &
                             grid(i)%eigvec_0,                          &
                             grid(i)%eigv,                              &
                             grid(i)%eigvec)
               endif
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
               WRITE(output,6) i
               WRITE(output,5)
         END IF
!
!              We are done with the DVR global construction.  Now we take
!              that information and construct the element matrices.
!
               DO j=1,num_reg(i)
                  n_alloc=nfun_reg(j,i)
                  WRITE(output,7) j, n_alloc
                  ALLOCATE(mat_reg_d(j,i)%pt_d(n_alloc),                  &
                           mat_reg_d(j,i)%nrm_d(n_alloc),                 &
                           mat_reg_d(j,i)%ke_mat_d(n_alloc,n_alloc),      &
                           mat_reg_d(j,i)%eigvec_mat_d(n_alloc,n_alloc),  &
                           mat_reg_d(j,i)%eigval_mat_d(n_alloc))
                  maxmem_reg=max(maxmem_reg,n_alloc)
               END DO        
               ALLOCATE(scr_d(5*maxmem_reg))
               IF(typke == 'dvr'.OR.typke == 'packed') THEN
                  mat_typ='full'
                  row_dim=nphy(i)
               ELSE
                  mat_typ='banded'
                  row_dim=row(i)
               END IF
!
!              This is a way to modify the element matrices to incorporate parts
!              of the diagonal which are time independent.
!
               IF(diag_mod /= 'none') THEN
                  call modify_diag(grid(i)%ke,grid(i)%v,row_dim,nphy(i),   &
                                   diag_mod,mat_typ)
               END IF
               start = 1
               DO j=1,num_reg(i)
                    call pt_reg(grid(i)%pt(start),                         &
                                grid(i)%wt(start),                         &
                                mat_reg_d(j,i)%pt_d,                       &
                                mat_reg_d(j,i)%nrm_d,                      &
                                nfun_reg(j,i),nphy(i),j)
                    start = start + nfun_reg(j,i) - 1
               END DO
!
!              Kinetic energy element matrices for this dimension
!
               start = 1
               IF(mat_typ=='full') THEN
                  DO j=1,num_reg(i)
                     call ke_reg_dvr(grid(i)%ke(start,start),              &
                                     mat_reg_d(j,i)%ke_mat_d,              &
                                     nfun_reg(j,i),nphy(i),j)
                     start = start + nfun_reg(j,i) - 1
                  END DO
               ELSE IF(mat_typ=='banded') THEN
                       DO j=1,num_reg(i)
                          call ke_reg_fd(grid(i)%ke(1,start),              &
                                         mat_reg_d(j,i)%ke_mat_d,          &
                                         nfun_reg(j,i),row_dim,j)
                          start = start + nfun_reg(j,i) - 1
                       END DO
               ELSE
                      CALL lnkerr('error')
               END IF
!
!              Diagonalize the KE matrices, one by one to get sector eigenvalues
!              and eigenvectors.
!
               DO j=1,num_reg(i)
!
!                 Diagonalize
!
                  call diag_reg(mat_reg_d(j,i)%ke_mat_d,                   &
                                mat_reg_d(j,i)%eigval_mat_d,               &
                                mat_reg_d(j,i)%eigvec_mat_d,               &
                                mat_reg_d(j,i)%eigvec_mat_d,               &
                                nfun_reg(j,i),j)
               END DO
               DEALLOCATE(scr_d)
      END DO
!
!              I have done NO DEALLOCATION as a user may want to link
!              directly to these matrices.  Since the code is terminated
!              if used as a standalone, they die anyway.
  ELSE
      WRITE(output,3)
      stop
  END IF
  stop
1    FORMAT(/,20X,'One-Body DVR Basis Code. coord = ',a24)
2    FORMAT(/,5x,'Nphy = ',i5,2x,'Nglobal = ',i5)
3    FORMAT(/,10x,'Error In Input File')
4    FORMAT(/15x,'Regional Information for DVR Basis,'                     &
                 ' Variable = ',i2)
5    FORMAT(/,10x,'Region',5x,'Number of Functions')            
6    FORMAT(/15x,'Regional Information for FD Basis,'             &
                 ' Variable = ',i2)
7    FORMAT(11x,i4,14x,i5)

END PROGRAM driver
