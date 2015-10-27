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
!
!
!
!     Compute the DVR points, weights, functions, 
!     first and second derivatives, kinetic energy matrix, 
!     eigenvalues and eigenvectors of the kinetic energy 
!     matrix, full one-particle Hamiltonian
!     matrix, eigenvalues and eigenvectors of the 
!     Hamiltonian matrix, one-body potential, 
!     value of DVR functions at the endpoints, 
!     and value of eigenvectors of the kinetic energy 
!     and one-body Hamiltonian at the endpoints.
!
  call dvr_basis(pt_0(1),                                   &
                 pt_n(1),                                   &
                 grid(1)%pt,                                &
                 grid(1)%wt,                                &
                 grid(1)%f,                                 &
                 grid(1)%df,                                &
                 grid(1)%ddf,                               &
                 grid(1)%ke,                                &
                 grid(1)%p_mom,                             &
                 grid(1)%eigv_0,                            &
                 grid(1)%eigvec_0,                          &
                 grid(1)%h,                                 &
                 grid(1)%eigv,                              &
                 grid(1)%eigvec,                            &
                 grid(1)%v,                                 &
                 grid(1)%srf_prm,                           &
                 grid(1)%srf_0,                             &
                 grid(1)%srf,                               &
                 coord(1),nphy(1),nglobal(1))
!
!
END subroutine Space_DVR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
