!deck Space_DVR.f
!***begin prologue     Space_DVR
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           space, propagation, dvr, finite difference
!***author             schneider, b. i.(nsf)
!***source
!***purpose            allocate and compute dvr functions and matrices.
!***
!***references
!***routines called    Name             Location
!                      ----             --------
!                    dvr_input        dvr_lib_f90
!                    dvr_basis        dvr_lib_f90
!                    fourier          dvr_lib_f90
!                    fourier_basis    dvr_lib_f90  
!                    fd_input         fd_lib_f90
!                    fd_basis         fd_lib_f90
!
!***modules used     Name               Location
!                    ----               --------
!                 dvrprop_global     dvr_lib_f90
!                 dvr_shared         dvr_lib_f90
!                 dvr_global         dvr_lib_f90
!
!***end prologue       space_prop  
  Subroutine Space_DVR
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  CHARACTER (LEN=1)                        :: itoc
  INTEGER                                  :: i, j
  INTEGER                                  :: start
!
! Allocate storage for the element matrices and potential
!   
  ALLOCATE(grid(spdim), num_reg(spdim), nfun_reg(maxreg,spdim) )
!
  DO  i=1,spdim
!
!-----------------------------------------------------------------------
!
!                   Begin the Code for the DVR Case
!
!-----------------------------------------------------------------------
!
!        Input the information needed to construct the DVR basis and matrix
!        elements for the i^{th} dimension.
!        At this point there is no explicit recognition of the 
!        "element" structure.  This will be done later.
!
      CALL dvr_input(nphy(i),nglobal(i),coord(i))
!
!         Allocate the needed DVR arrays.
!
      IF(diag) then                      
         ALLOCATE(grid(i)%eigv_0(nphy(i)),                         &   
                  grid(i)%eigvec_0(nphy(i),nphy(i)),               &   
                  grid(i)%eigv(nphy(i)),                           &   
                  grid(i)%eigvec(nphy(i),nphy(i)),                 &   
                  grid(i)%srf_0(nphy(i),2),                        &
                  grid(i)%srf(nphy(i),2))    
      END IF
!
      IF(typwt /= 'fourier') THEN
         ALLOCATE(grid(i)%pt(nphy(i)),                             &
                  grid(i)%wt(nphy(i)),                             &
                  grid(i)%f(nphy(i),nphy(i)),                      &
                  grid(i)%df(nphy(i),nphy(i)),                     &   
                  grid(i)%ddf(nphy(i),nphy(i)),                    &   
                  grid(i)%ke(nphy(i),nphy(i)),                     &   
                  grid(i)%p_mom(nphy(i),nphy(i)),                  &   
                  grid(i)%h(nphy(i),nphy(i)),                      &   
                  grid(i)%v(nphy(i)),grid(i)%srf_prm(2))
!
!           Compute the DVR points, weights, functions, 
!           first and second derivatives, kinetic energy matrix, 
!           eigenvalues and eigenvectors of the kinetic energy 
!           matrix, full one-particle Hamiltonian
!           matrix, eigenvalues and eigenvectors of the 
!           Hamiltonian matrix, one-body potential, 
!           value of DVR functions at the endpoints, 
!           and value of eigenvectors of the kinetic energy 
!           and one-body Hamiltonian at the endpoints.
!
         call dvr_basis(pt_0(i),                                   &
                        pt_n(i),                                   &
                        grid(i)%pt,                                &
                        grid(i)%wt,                                &
                        grid(i)%f,                                 &
                        grid(i)%df,                                &
                        grid(i)%ddf,                               &
                        grid(i)%ke,                                &
                        grid(i)%p_mom,                             &
                        grid(i)%eigv_0,                            &
                        grid(i)%eigvec_0,                          &
                        grid(i)%h,                                 &
                        grid(i)%eigv,                              &
                        grid(i)%eigvec,                            &
                        grid(i)%v,                                 &
                        grid(i)%srf_prm,                           &
                        grid(i)%srf_0,                             &
                        grid(i)%srf,                               &
                        coord(i),nphy(i),nglobal(i))
!
!        Deallocate everything except the points, weights, 
!        functions, kinetic energy and potential arrays
!
         DEALLOCATE(grid(i)%df,grid(i)%ddf,                        &   
                    grid(i)%p_mom,grid(i)%h,grid(i)%srf_prm)
      ELSE IF(typwt == 'fourier') THEN
         ALLOCATE(grid(i)%pt(nphy(i)),                             &
                  grid(i)%wt(nphy(i)),                             &
                  grid(i)%f(nphy(i),nphy(i)),                      &
                  grid(i)%ke(nphy(i),nphy(i)),                     &   
                  grid(i)%h(nphy(i),nphy(i)),                      &   
                  grid(i)%v(nphy(i)))
         CALL fourier(grid(i)%pt,grid(i)%wt,grid(i)%f,             &
                      edge,nphy(i))
         call fourier_basis(grid(i)%pt,                            &
                            grid(i)%ke,                            &
                            grid(i)%h,                             &
                            grid(i)%v,                             &
                            grid(i)%eigv_0,                        &
                            grid(i)%eigvec_0,                      &
                            grid(i)%eigv,                          &
                            grid(i)%eigvec,                        &
                            coord(i),nphy(i))
!
         DEALLOCATE(grid(i)%h)
      ELSE
         call lnkerr('weight error')
      END IF
      row(i) = 2*nphy(i) - 2
      IF(diag) then                      
         DEALLOCATE(grid(i)%srf_0,grid(i)%srf)
      END IF
      num_reg(i) = nreg
      IF(bcl == 0) then
         npt(1) = npt(1) - 1
      END IF
      IF(bcr == 0 ) then
         npt(nreg) = npt(nreg) - 1
      END IF
      nfun_reg(:,i) = npt
  END DO
  n_reg_real = num_reg 
  n3d=1
  DO  i=1,spdim
      n3d=n3d*nphy(i)
  END DO
  key='FEDVR'
  WRITE(iout,1) spdim, key, n3d
!
!
1 Format(/,1x,'Number of Dimensions   = ',i2,   &
           1x,'Spatial discretization = ',a8,   &
         /,1x,'Number of Points       = ',i10)  
END subroutine Space_DVR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
