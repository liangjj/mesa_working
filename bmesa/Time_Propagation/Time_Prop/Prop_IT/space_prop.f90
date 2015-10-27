!deck space_prop.f
!***begin prologue     space_prop
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           space, propagation, dvr, finite difference
!***author             schneider, b. i.(nsf)
!***source
!***purpose            allocate and compute dvr functions and matrices.
!***
!***references
!***routines called    iosys, util, mdutil, dvr_input, dvr_basis, fourier
!***                   fourier_basis, fd_input, fd_basis,
!***modules            dvrprop_global_it, dvr_shared, dvr_global
!***end prologue       space_prop  
  Subroutine Space_Prop
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
      IF(typke == 'dvr'.OR.typke == 'packed') THEN
         CALL dvr_input(nphy(i),nglobal(i),coord(i))
!
!         Allocate the needed DVR arrays.
!
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
            call iosys('write real "eigenvalues for variable-'        &
                                    //itoc(i)//'" to bec',nphy(i),    &
                                      grid(i)%eigv,0,' ')
            call iosys('write real "eigenvectors for variable-'       &
                                    //itoc(i)//'" to bec',            &
                                      nphy(i)*nphy(i),                &
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
         ALLOCATE(grid(i)%pt(nphy(i)),                                &
                  grid(i)%wt(nphy(i)),                                &
                  grid(i)%ke(row(i),nphy(i)),                         &
                  grid(i)%h(row(i),nphy(i)),                          &
                  grid(i)%v(nphy(i)))
         IF(diag) then
            ALLOCATE(grid(i)%eigv_0(nphy(i)),                         &
                     grid(i)%eigvec_0(nphy(i),nphy(i)),               &
                     grid(i)%eigv(nphy(i)),                           &
                     grid(i)%eigvec(nphy(i),nphy(i)))
         END IF
         CALL fd_basis(pt_0(i),                                       &
                       grid(i)%pt,                                    &
                       grid(i)%wt,                                    &
                       grid(i)%ke,                                    &
                       grid(i)%eigv_0,                                &
                       grid(i)%eigvec_0,                              &
                       grid(i)%h,                                     &
                       grid(i)%eigv,                                  &
                       grid(i)%eigvec,                                &
                       grid(i)%v,                                     &
                       nphy(i),nglobal(i),row(i),coord(i))
         DEALLOCATE(grid(i)%h)
         IF(diag) THEN
!
!           Save one body eigenvalues and vectors.
!
          
            call iosys('write real "eigenvalues for variable-'        &
                                    //itoc(i)//'" to bec',nphy(i),    &
                                      grid(i)%eigv,0,' ')
            call iosys('write real "eigenvectors for variable-'       &
                                    //itoc(i)//'" to bec',            &
                                      nphy(i)*nphy(i),                &
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
      END IF
  END DO
  n3d=1
  DO  i=1,spdim
      n3d=n3d*nphy(i)
  END DO
END subroutine Space_Prop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
