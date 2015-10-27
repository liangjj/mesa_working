!deck CC_Space_Prop.f
!***begin prologue     CC_Space_Prop
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
!
!***modules used       Name             Location
!                      ----             --------
!                   dvrprop_global     dvr_lib_f90
!                   dvr_shared         dvr_lib_f90
!                   dvr_global         dvr_lib_f90
!
!***end prologue       CC_Space_Prop  
  Subroutine CC_Space_Prop
  USE dvrprop_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
!
! Allocate storage for the element matrices and potential
! There is one grid for all the channels
!  
  ALLOCATE(grid, num_reg, nfun_reg(maxreg) )
!
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

  CALL dvr_input(nphy,nglobal,coord)
!
!      Allocate the needed DVR arrays.
!
!
  IF(diag) then                      
     ALLOCATE(grid%eigv_0(nphy),                      &   
              grid%eigvec_0(nphy,nphy),               &   
              grid%eigv(nphy),                        &   
              grid%eigvec(nphy,nphy),                 &   
              grid%srf_0(nphy,2),                     &
              grid%srf(nphy,2))    
  END IF
!
  ALLOCATE(grid%pt(nphy),                             &
           grid%wt(nphy),                             &
           grid%f(nphy,nphy),                         &
           grid%df(nphy,nphy),                        &   
           grid%ddf(nphy,nphy),                       &   
           grid%ke(nphy,nphy),                        &   
           grid%p_mom(nphy,nphy),                     &   
           grid%h(nphy,nphy),                         &   
           grid%v(nphy),grid%srf_prm(2))
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
  CALL dvr_basis(pt_0,                                &
                 grid%pt,                             &
                 grid%wt,                             &
                 grid%f,                              &
                 grid%df,                             &
                 grid%ddf,                            &
                 grid%ke,                             &
                 grid%p_mom,                          &
                 grid%eigv_0,                         &
                 grid%eigvec_0,                       &
                 grid%h,                              &
                 grid%eigv,                           &
                 grid%eigvec,                         &
                 grid%v,                              &
                 grid%srf_prm,                        &
                 grid%srf_0,                          &
                 grid%srf,                            &
                 coord,nphy,nglobal)
!
!        Deallocate everything except the points, weights, 
!        functions, kinetic energy and potential arrays
!
  DEALLOCATE(grid%df,grid%ddf,                        &   
             grid%p_mom,grid%h,grid%srf_prm)
  row = 2*nphy - 2
  IF(diag) then                      
     DEALLOCATE(grid%srf_0,grid%srf)
  END IF
  num_reg = nreg
  IF(bcl == 0) then
     npt(1) = npt(1) - 1
  END IF
  IF(bcr == 0 ) then
     npt(nreg) = npt(nreg) - 1
  END IF
  nfun_reg = npt
  n_reg_real = num_reg -  n_reg_absorb  
  vec_len=nchan*nphy
  key='FEDVR'
  WRITE(iout,1) vec_len
!
!
1 Format(/,1x,'Vector Length = ',i10)  
END subroutine CC_Space_Prop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
