!***********************************************************************
! regional_module
!**begin prologue     regional_module
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains all of the subroutines to form,
!***                  diagonalize and construct the regional matrices.
!***                  Explicit interfaces are used to allow
!***                  a transparent use of generic subroutines which work
!***                  for both real and complex vectors.  This feature
!***                  permits a single code to be used for both real and
!***                  imaginary time propagation.
!***description       See subroutines
!**references
!**modules needed     See USE statements below
!**end prologue       regional_module
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      MODULE regional_module
                      USE dvr_global
                      USE dvr_shared
                      USE dvrprop_global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck regional_matrices.f
!***begin prologue     regional_matrices
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           space, propagation, dvr, finite difference
!***author             schneider, b. i.(nsf)
!***source
!***purpose            the global FEDVR matrices are split into their
!***                   regional parts for ease of use in the propagation.
!***                   the regional kinetic energy matrix may be modified
!***                   to include any diagonal contribution from the one
!***                   body potential or any other diagnal part.  once this
!***                   is done, the kinetic energy is diagonalized.  the 
!***                   eigenvectors and eigenvalues are used to form the
!***                   element propagators later in the code.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       space
  Subroutine Regional_matrices
  IMPLICIT NONE
  INTEGER                                  :: i, j, start, row_dim
  CHARACTER(LEN=8)                         :: mat_typ
!
!
  ALLOCATE(mat_reg(maxreg,spdim))
  IF(algorithm == 'split_operator') THEN
     DO  i=1,spdim
         IF(typke == 'dvr'.OR.typke == 'packed') THEN
            write(iout,1) i
         ELSE
            write(iout,2) i
         END IF
         WRITE(iout,3)
         maxmem_reg=0
         DO j=1,n_reg_real(i)
            write(iout,4) j, nfun_reg(j,i), 'real'
            ALLOCATE(mat_reg(j,i)%pt_d                             &
                    (nfun_reg(j,i)),                               &
                     mat_reg(j,i)%ke_mat_d                         &
                    (nfun_reg(j,i),nfun_reg(j,i)),                 &
                     mat_reg(j,i)%eigval_mat_d                     &
                    (nfun_reg(j,i)),                               &
                     mat_reg(j,i)%eigvec_mat_d                     &
                    (nfun_reg(j,i),nfun_reg(j,i)))
            maxmem_reg=max(maxmem_reg,nfun_reg(j,i))
         END DO
         ALLOCATE(scr_d(5*maxmem_reg))
         IF(absorb) THEN
!
!           There is an absorber.  Allocate the complex arrays.
!
            DO j=n_reg_real(i)+1, num_reg(i)
               write(iout,4) j, nfun_reg(j,i),'complex'
               ALLOCATE(mat_reg(j,i)%pt_d                          &
                       (nfun_reg(j,i)),                            &
                        mat_reg(j,i)%v_add_z                       &
                       (nfun_reg(j,i)),                            &
                        mat_reg(j,i)%ke_mat_d                      &
                       (nfun_reg(j,i),nfun_reg(j,i)),              &
                        mat_reg(j,i)%ke_mat_z                      &
                       (nfun_reg(j,i),nfun_reg(j,i)),              &
                        mat_reg(j,i)%eigval_mat_z                  &
                       (nfun_reg(j,i)),                            &
                        mat_reg(j,i)%eigvec_mat_z_r                &
                       (nfun_reg(j,i),nfun_reg(j,i)),              &
                        mat_reg(j,i)%eigvec_mat_z_l                &
                       (nfun_reg(j,i),nfun_reg(j,i)))
               maxmem_reg=max(maxmem_reg,nfun_reg(j,i))
            END DO
            ALLOCATE(scr_z(5*maxmem_reg))
         END IF
         IF(typke == 'dvr'.OR.typke == 'packed') THEN
            mat_typ='full'
            row_dim=nphy(i)
         ELSE
            mat_typ='banded'
            row_dim=row(i)
         END IF
         IF(diag_mod /= 'none') THEN
            call modify_diag(grid(i)%ke,grid(i)%v,row_dim,nphy(i), &
                             diag_mod,mat_typ)
         END IF
         start = 1
         DO j=1,num_reg(i)
            call pt_reg(grid(i)%pt(start),mat_reg(j,i)%pt_d,       &
                        nfun_reg(j,i),nphy(i),j)
            start = start + nfun_reg(j,i) - 1
         END DO
!
!        here we actually form the FEDVR kinetic energy regional
!        matrices.
!
         start = 1
         IF(mat_typ=='full') THEN
            DO j=1,num_reg(i)
               call ke_reg_dvr(grid(i)%ke(start,start),            &
                               mat_reg(j,i)%ke_mat_d,              &
                               nfun_reg(j,i),nphy(i),j)
               start = start + nfun_reg(j,i) - 1
            END DO
         ELSE IF(mat_typ=='banded') THEN
            DO j=1,num_reg(i)
               call ke_reg_fd(grid(i)%ke(1,start),                 &
                              mat_reg(j,i)%ke_mat_d,               &
                              nfun_reg(j,i),row_dim,j)
               start = start + nfun_reg(j,i) - 1
            END DO
         ELSE
            CALL lnkerr('error')
         END IF  
         IF(absorb) THEN
            DO j=n_reg_real(i)+1,num_reg(i)
               call add_absorb(mat_reg(j,i)%pt_d,                  &
                               mat_reg(j,i)%ke_mat_d,              &
                               mat_reg(j,i)%v_add_z,               &
                               mat_reg(j,i)%ke_mat_z,              &
                               nfun_reg(j,i))
            END DO
         END IF
         DO j=1,n_reg_real(i)
!
!           Diagonalize
!
            call regional_diag_d(mat_reg(j,i)%ke_mat_d,            &
                                 mat_reg(j,i)%eigval_mat_d,        &
                                 mat_reg(j,i)%eigvec_mat_d,        &
                                 nfun_reg(j,i),j)
         END DO
         DEALLOCATE(scr_d)
         IF(absorb) THEN
            DO j=n_reg_real(i)+1,num_reg(i)
!
!           Diagonalize
!
               call regional_diag_z(mat_reg(j,i)%ke_mat_z,         &
                                    mat_reg(j,i)%eigval_mat_z,     &
                                    mat_reg(j,i)%eigvec_mat_z_r,   &
                                    mat_reg(j,i)%eigvec_mat_z_l,   &
                                    nfun_reg(j,i),j)
            END DO
            DEALLOCATE(scr_z)
         END IF
     END DO
  ELSE
     DO  i=1,spdim
         IF(typke == 'dvr'.OR.typke == 'packed') THEN
            write(iout,1) i
         ELSE
            write(iout,2) i
         END IF
         WRITE(iout,3)
         DO j=1,n_reg_real(i)
            write(iout,4) j, nfun_reg(j,i), 'real'
            ALLOCATE(mat_reg(j,i)%pt_d                             &
                    (nfun_reg(j,i)),                               &
                     mat_reg(j,i)%ke_mat_d                         &
                    (nfun_reg(j,i),nfun_reg(j,i)))
         END DO
         IF(absorb) THEN
!
!           There is an absorber.  Allocate the complex arrays.
!
            DO j=n_reg_real(i)+1, num_reg(i)
               write(iout,4) j, nfun_reg(j,i),'complex'
               ALLOCATE(mat_reg(j,i)%pt_d                          &
                       (nfun_reg(j,i)),                            &
                        mat_reg(j,i)%v_add_z                       &
                       (nfun_reg(j,i)),                            &
                        mat_reg(j,i)%ke_mat_d                      &
                       (nfun_reg(j,i),nfun_reg(j,i)),              &
                        mat_reg(j,i)%ke_mat_z                      &
                       (nfun_reg(j,i),nfun_reg(j,i)) )
            END DO
         END IF
         IF(typke == 'dvr'.OR.typke == 'packed') THEN
            mat_typ='full'
            row_dim=nphy(i)
         ELSE
            mat_typ='banded'
            row_dim=row(i)
        END IF
       IF(diag_mod /= 'none') THEN
            call modify_diag(grid(i)%ke,grid(i)%v,row_dim,nphy(i), &
                             diag_mod,mat_typ)
       END IF
       start = 1
       DO j=1,num_reg(i)
          call pt_reg(grid(i)%pt(start),mat_reg(j,i)%pt_d,       &
                      nfun_reg(j,i),nphy(i),j)
          start = start + nfun_reg(j,i) - 1
       END DO
       start = 1
       IF(mat_typ=='full') THEN
          DO j=1,num_reg(i)
             call ke_reg_dvr(grid(i)%ke(start,start),            &
                             mat_reg(j,i)%ke_mat_d,              &
                             nfun_reg(j,i),nphy(i),j)
             start = start + nfun_reg(j,i) - 1
          END DO
       ELSE IF(mat_typ=='banded') THEN
          DO j=1,num_reg(i)
             call ke_reg_fd(grid(i)%ke(1,start),                 &
                            mat_reg(j,i)%ke_mat_d,               &
                            nfun_reg(j,i),row_dim,j)
             start = start + nfun_reg(j,i) - 1
          END DO
       ELSE
          CALL lnkerr('error')
       END IF  
       IF(absorb) THEN
          DO j=n_reg_real(i)+1,num_reg(i)
             call add_absorb(mat_reg(j,i)%pt_d,                  &
                             mat_reg(j,i)%ke_mat_d,              &
                             mat_reg(j,i)%v_add_z,               &
                             mat_reg(j,i)%ke_mat_z,              &
                             nfun_reg(j,i))
          END DO
       END IF
       IF(absorb) THEN
          DO j=n_reg_real(i) + 1, num_reg(i) 
          END DO
       END IF
     END DO
  END IF
1 FORMAT(/15x,'Regional Information for DVR Basis,'                &
              ' Variable = ',i2)
2 FORMAT(/15x,'Regional Information for FD Basis,'                 &
                 ' Variable = ',i2)
3 FORMAT(/,17x,'Region',5x,'Number of Functions',5x,'Type')
4 FORMAT(16x,i4,14x,i5,9x,a8)
END Subroutine regional_matrices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck regional_diag_d
!***begin prologue     regional_diag_d
!***date written       040706   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Regional matrix elements of kinetic energy
!                      with zero diagonal are diagonalized.
!***
!***references
!***routines called    dsyev(LAPACK)
!***end prologue       regional_diag_d
!
  SUBROUTINE regional_diag_d(k_region,eval_region,evec_region,         &
                           n_region,reg)
  IMPLICIT NONE
  INTEGER                                :: n_region, reg
  INTEGER                                :: info
  REAL*8, DIMENSION(n_region,n_region)   :: k_region
  REAL*8, DIMENSION(n_region)            :: eval_region
  REAL*8, DIMENSION(n_region,n_region)   :: evec_region
  CHARACTER(LEN=4)                       :: itoc
!
!
  evec_region = k_region
  call dsyev('v','l',n_region,evec_region,n_region,eval_region,   &
              scr_d,5*n_region,info)
  if(log_main(1)) then
     title='sector eigenvalues region = '//itoc(reg)
     call prntfmn(title,eval_region,n_region,1,n_region,1,iout,'e')
     title='sector eigenvectors = '//itoc(reg)
     call prntfmn(title,evec_region,n_region,n_region,            &
                                    n_region,n_region,iout,'e')
  end if
  END SUBROUTINE regional_diag_d  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck regional_diag_z
!***begin prologue     regional_diag_z
!***date written       040706   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            regional matrix elements of kinetic energy
!                      with complex potential is diagonalized.
!***
!***references
!***routines called    dsyev(LAPACK)
!***end prologue       regional_diag_z
!
  SUBROUTINE regional_diag_z(k_region,eval_region,evec_region_right,    &
                             evec_region_left,n_region,reg)
  IMPLICIT NONE
  INTEGER                                  :: n_region, reg
  INTEGER                                  :: info
  INTEGER                                  :: i, j
  COMPLEX*16, DIMENSION(n_region,n_region) :: k_region
  COMPLEX*16, DIMENSION(n_region)          :: eval_region
  COMPLEX*16, DIMENSION(n_region,n_region) :: evec_region_right
  COMPLEX*16, DIMENSION(n_region,n_region) :: evec_region_left
  COMPLEX*16                               :: val_z, cdotc
!
!
  call zgeev('v','v',n_region,k_region,n_region,eval_region,       &
             evec_region_left,n_region,evec_region_right,n_region, &
             scr_z,10*n_region,scr_d,info)
  IF(info /= 0 ) THEN
     CALL lnkerr('eigenvalue routine error')
  END IF
!
! Renormalize
!
  DO i=1,n_region
     val_z=1.d0/sqrt(cdotc(n_region,evec_region_left(1,i),1,       &
                     evec_region_right(1,i),1))   
     DO j=1,n_region
        evec_region_right(j,i) = val_z * evec_region_right(j,i)
        evec_region_left(j,i) = conjg(val_z) * evec_region_left(j,i)
     END DO
  END DO
  title='sector eigenvalues complex region'
  call prntcmn(title,eval_region,n_region,1,n_region,1,iout,'e')
  if(log_main(1)) then
     title='sector left eigenvectors'
     call prntcmn(title,evec_region_left,n_region,n_region,        &
                  n_region,n_region,iout,'e')
     title='sector right eigenvectors'
     call prntcmn(title,evec_region_right,n_region,n_region,       &
                  n_region,n_region,iout,'e')
  end if
END SUBROUTINE regional_diag_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck regional_prop_d.f
!***begin prologue     regional_prop_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            propagators
!***
!***description        The second order decomposition of the time exponential
!***                   is performed for a Hamiltonian H = H_a + H_b + V as
!***                   U_2(tau) = exp(-V * tau/2)   *
!***                                   u_2(tau)     *
!***                              exp(-V * tau/2)
!***
!***                   u_2(tau) = exp(-H_a * tau/2) * 
!***                              exp(-H_b * tau)   * 
!***                              exp(_H_a * tau/2) *
!***
!***                   The expression for u_2 requires the calculation of two 
!***                   independent propagators.
!***
!***                   The fourth order decomposition of the time exponential
!***                   is performed as
!***                   U_4(tau) = U_2(p*tau)       * 
!***                              U_2(p*tau)       * 
!***                              U_2((1-4*p)*tau) * 
!***                              U_2(p*tau)       * 
!***                              U_2(p*tau) 
!***                              p = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
!***
!***                   The expression for U_4 requires the calculation of four 
!***                   independent propagators.
!***references
!***routines called
!***end prologue       propagator
!
  SUBROUTINE regional_prop_d(q,t_in)
  IMPLICIT NONE
  INTEGER                                  :: q, i
  INTEGER                                  :: trips, ntrip
  REAL*8                                   :: t_in, tau
!
!
! Special case of only one region
!
  write(iout,1) prop_order, p_fac
!
  IF ( n_reg == 1) THEN
       write(iout,2)
       IF(prop_order==2) THEN
          tau = t_in
          call umat_reg_d(mat_reg(starting_reg,q)%eigval_mat_d,                    &
                          mat_reg(starting_reg,q)%eigvec_mat_d,                    &
                          mat_reg(starting_reg,q)%exp_t_mat_d(:,:,1),              &
                          exp_d,tau,nfun_reg(starting_reg,q))
       ELSE IF(prop_order==4) THEN
          p_loc = p_fac
          trips = 2
          Write(iout,3) trips
          DO ntrip=1,trips
             write(iout,4) ntrip, p_loc
              tau = p_loc * t_in
              call umat_reg_d(mat_reg(starting_reg,q)%eigval_mat_d,                &
                              mat_reg(starting_reg,q)%eigvec_mat_d,                &
                              mat_reg(starting_reg,q)%exp_t_mat_d(:,:,ntrip),      &
                              exp_d,tau,nfun_reg(starting_reg,q))
!                                                                          
!             Modify p_loc for fourth order propagator so it is 
!             correct second time through last loop.            
!                                                                          
              p_loc = ( 1.d0 - 4.d0 * p_loc )
          END DO
       END IF
  ELSE
       p_loc = p_fac
       trips = 1
       IF ( prop_order == 4 ) THEN
            trips = 2
       END IF
       Write(iout,3) trips
!
       DO ntrip=1,trips
          write(iout,4) ntrip, p_loc
!
!         Construct the "outer" propagator at p_fac * delta_t/2
!
          tau = p_loc * t_in *.5d0
          DO i = starting_reg, ending_reg , 2
             WRITE(iout,5) i
             call umat_reg_d(mat_reg(i,q)%eigval_mat_d,                            &
                             mat_reg(i,q)%eigvec_mat_d,                            &
                             mat_reg(i,q)%exp_t_mat_d(:,:,ntrip),                  &
                             exp_d,tau,nfun_reg(i,q))
          END DO
!
!         Construct the "inner" propagator at p_loc * delta_t
!
          tau = p_loc * t_in
          DO i=starting_reg + 1 , ending_reg , 2
             WRITE(iout,5) i
             call umat_reg_d(mat_reg(i,q)%eigval_mat_d,                            &
                             mat_reg(i,q)%eigvec_mat_d,                            &
                             mat_reg(i,q)%exp_t_mat_d(:,:,ntrip),                  &
                             exp_d,tau,nfun_reg(i,q))
          END DO
!
!         Modify p_loc for fourth order propagator so it's correct second time 
!         through last loop.
!
          p_loc = ( 1.d0 - 4.d0 * p_loc )
       END DO
  END IF
1 Format(/,10x,'Constructing the Propagators for order = ',i2, &
         /,10x,'Exponential factor                     = ',e15.8)
2 FORMAT(/,10x,'There is only one region')
3 FORMAT(/,10x,'Number of Propagator Passes            = ',i1)
4 FORMAT(/,10x,'Pass                                   = ',i1, &
         /,10x,'p_fac                                  = ',e15.8)
5 FORMAT(/,5x,'Sector = ',i5)
END SUBROUTINE regional_prop_d
!***********************************************************************
!***********************************************************************
!deck regional_prop_h.f
!***begin prologue     regional_prop_h
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            propagators
!***
!***description        The second order decomposition of the time exponential
!***                   is performed for a Hamiltonian H = H_a + H_b + V as
!***                   U_2(tau) = exp(-V * tau/2)   *
!***                                   u_2(tau)     *
!***                              exp(-V * tau/2)
!***
!***                   u_2(tau) = exp(-H_a * tau/2) * 
!***                              exp(-H_b * tau)   * 
!***                              exp(_H_a * tau/2) *
!***
!***                   The expression for u_2 requires the calculation of two 
!***                   independent propagators.
!***
!***                   The fourth order decomposition of the time exponential
!***                   is performed as
!***                   U_4(tau) = U_2(p*tau)       * 
!***                              U_2(p*tau)       * 
!***                              U_2((1-4*p)*tau) * 
!***                              U_2(p*tau)       * 
!***                              U_2(p*tau) 
!***                              p = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
!***
!***                   The expression for U_4 requires the calculation of four 
!***                   independent propagators.
!***references
!***routines called
!***end prologue       propagator
!
  SUBROUTINE regional_prop_h(q,t_in)
  IMPLICIT NONE
  INTEGER                                  :: q, i
  INTEGER                                  :: trips, ntrip
  REAL*8                                   :: t_in, tau
!
!
! Special case of only one region
!
  write(iout,1) prop_order, p_fac
!
  IF ( n_reg == 1) THEN
       write(iout,2)
       IF(prop_order==2) THEN
          tau = t_in
          call umat_reg_h(mat_reg(starting_reg,q)%eigval_mat_d,                    &
                          mat_reg(starting_reg,q)%eigvec_mat_d,                    &
                          mat_reg(starting_reg,q)%exp_t_mat_z(:,:,1),              &
                          exp_z,tau,nfun_reg(starting_reg,q))
       ELSE IF(prop_order==4) THEN
          p_loc = p_fac
          trips = 2
          Write(iout,3) trips
          DO ntrip=1,trips
             write(iout,4) ntrip, p_loc
              tau = p_loc * t_in
              call umat_reg_h(mat_reg(starting_reg,q)%eigval_mat_d,                &
                              mat_reg(starting_reg,q)%eigvec_mat_d,                &
                              mat_reg(starting_reg,q)%exp_t_mat_z(:,:,ntrip),      &
                              exp_z,tau,nfun_reg(starting_reg,q))
!                                                                          
!             Modify p_loc for fourth order propagator so it is 
!             correct second time through last loop.            
!                                                                          
              p_loc = ( 1.d0 - 4.d0 * p_loc )
          END DO
       END IF
  ELSE
       p_loc = p_fac
       trips = 1
       IF ( prop_order == 4 ) THEN
            trips = 2
       END IF
       Write(iout,3) trips
!
       DO ntrip=1,trips
          write(iout,4) ntrip, p_loc
!
!         Construct the "outer" propagator at p_fac * delta_t/2
!
          tau = p_loc * t_in *.5d0
          DO i = starting_reg, ending_reg , 2
             WRITE(iout,5) i
             call umat_reg_h(mat_reg(i,q)%eigval_mat_d,                            &
                             mat_reg(i,q)%eigvec_mat_d,                            &
                             mat_reg(i,q)%exp_t_mat_z(:,:,ntrip),                  &
                             exp_z,tau,nfun_reg(i,q))
          END DO
!
!         Construct the "inner" propagator at p_loc * delta_t
!
          tau = p_loc * t_in
          DO i=starting_reg + 1 , ending_reg , 2
             WRITE(iout,5) i
             call umat_reg_h(mat_reg(i,q)%eigval_mat_d,                            &
                             mat_reg(i,q)%eigvec_mat_d,                            &
                             mat_reg(i,q)%exp_t_mat_z(:,:,ntrip),                  &
                             exp_z,tau,nfun_reg(i,q))
          END DO
!
!         Modify p_loc for fourth order propagator so it's correct second time 
!         through last loop.
!
          p_loc = ( 1.d0 - 4.d0 * p_loc )
       END DO
  END IF
1 Format(/,10x,'Constructing the Propagators for order = ',i2, &
         /,10x,'Exponential factor                     = ',e15.8)
2 FORMAT(/,10x,'There is only one region')
3 FORMAT(/,10x,'Number of Propagator Passes            = ',i1)
4 FORMAT(/,10x,'Pass                                   = ',i1, &
         /,10x,'p_fac                                  = ',e15.8)
5 FORMAT(/,5x,'Sector = ',i5)
END SUBROUTINE regional_prop_h
!
!***********************************************************************
!***********************************************************************
  SUBROUTINE regional_prop_z(q,t_in)
  IMPLICIT NONE
  INTEGER                                  :: q
  INTEGER                                  :: i
  INTEGER                                  :: trips, ntrip
  REAL*8                                   :: t_in, tau
!
!
  write(iout,1) prop_order, p_fac
!
!-----------------------------------------------------------------------
!
! Special case of only one region
!
  IF ( n_reg == 1) THEN
       write(iout,2)
       IF(prop_order==2) then
          tau = t_in
          call umat_reg_z(mat_reg(starting_reg,q)%eigval_mat_z,                    &
                          mat_reg(starting_reg,q)%eigvec_mat_z_r,                  &

                          mat_reg(starting_reg,q)%eigvec_mat_z_l,                  &
                          mat_reg(starting_reg,q)%exp_t_mat_z(:,:,1),              &
                          exp_z,tau,nfun_reg(starting_reg,q))
       ELSE IF(prop_order==4) THEN
          p_loc = p_fac
          trips = 2
          Write(iout,3) trips
          DO ntrip=1,trips
             write(iout,4) ntrip, p_loc
              tau = p_loc * t_in
              call umat_reg_z(mat_reg(starting_reg,q)%eigval_mat_z,                &
                              mat_reg(starting_reg,q)%eigvec_mat_z_r,              &
                              mat_reg(starting_reg,q)%eigvec_mat_z_l,              &
                              mat_reg(starting_reg,q)%exp_t_mat_z(:,:,ntrip),      &
                              exp_z,tau,nfun_reg(starting_reg,q))
!
!             Modify p_loc for fourth order propagator so it is correct second time
!             through last loop.
!
              p_loc = ( 1.d0 - 4.d0 * p_loc )
          END DO
       END IF
  ELSE
       p_loc = p_fac
       trips = 1
       IF ( prop_order == 4 ) THEN
            trips = 2
       END IF
       Write(iout,3) trips
!
       DO ntrip=1,trips
          write(iout,4) ntrip, p_loc
!
!         Construct the "outer" propagator at p_loc * delta_t/2
!
          tau = p_loc * t_in *.5d0
          DO i=starting_reg, ending_reg, 2
             call umat_reg_z(mat_reg(i,q)%eigval_mat_z,                            &
                             mat_reg(i,q)%eigvec_mat_z_r,                          &
                             mat_reg(i,q)%eigvec_mat_z_l,                          &
                             mat_reg(i,q)%exp_t_mat_z(:,:,ntrip),                  &
                             exp_z,tau,nfun_reg(i,q))
          END DO
!
!         Construct the "inner" propagator at p_loc * delta_t
!
          tau = p_loc * t_in
          DO i=starting_reg + 1, ending_reg, 2
             call umat_reg_z(mat_reg(i,q)%eigval_mat_z,                            &
                             mat_reg(i,q)%eigvec_mat_z_r,                          &
                             mat_reg(i,q)%eigvec_mat_z_l,                          &
                             mat_reg(i,q)%exp_t_mat_z(:,:,ntrip),                  &
                             exp_z,tau,nfun_reg(i,q))
          END DO
!
!         Modify p_loc for fourth order propagator so it is correct second time
!         through last loop.
!
          p_loc = ( 1.d0 - 4.d0 * p_loc )
       END DO
  END IF
1 Format(/,10x,'Constructing the Propagators for order = ',i2, &
         /,10x,'Exponential factor                     = ',e15.8)
2 FORMAT(/,10x,'There is only one region')
3 FORMAT(/,10x,'Number of Propagator Passes            = ',i1)
4 FORMAT(/,10x,'Pass                                   = ',i1, &
         /,10x,'p_fac                                  = ',e15.8)
5 FORMAT(/,5x,'Sector = ',i5)
END SUBROUTINE regional_prop_z
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**begin prologue     regional_umat
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       regional_umat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  SUBROUTINE umat_reg_d(eval_region,evec_region,exp_it,exp_it_d,tau,n_region)
  IMPLICIT NONE
  INTEGER                                     :: n_region
  INTEGER                                     :: i
  REAL*8, DIMENSION(n_region)                 :: eval_region
  REAL*8, DIMENSION(n_region,n_region)        :: evec_region
  REAL*8, DIMENSION(n_region,n_region)        :: exp_it
  REAL*8, DIMENSION(n_region,n_region)        :: exp_it_d
  REAL*8                                      :: ex_fac, tau
  CHARACTER*16                                :: fptoc
  ex_fac = tau/hbar
  DO i=1,n_region
     exp_it_d(:,i) = evec_region(:,i) * exp(- eval_region(i) * ex_fac)  
  END DO
  call ebct(exp_it,exp_it_d,evec_region,n_region,n_region,n_region)
  IF (log_main(1)) then
      title='exponential sector propagator at tau = '//fptoc(tau)
      call prntfmn(title,exp_it,n_region,n_region,                 &
                   n_region,n_region,iout,'e')
  END IF
  END SUBROUTINE umat_reg_d 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***begin prologue     umat_reg_h
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            regional unitary propagators.
!***
!***references
!***routines called
!***end prologue       umat_reg_h
!
  SUBROUTINE umat_reg_h(eval_region,evec_region,exp_rt,exp_rt_z,tau,n_region)
  IMPLICIT NONE
  INTEGER                                     :: n_region
  INTEGER                                     :: i
  REAL*8, DIMENSION(n_region)                 :: eval_region
  REAL*8, DIMENSION(n_region,n_region)        :: evec_region
  COMPLEX*16, DIMENSION(n_region,n_region)    :: exp_rt
  COMPLEX*16, DIMENSION(n_region,n_region)    :: exp_rt_z
  REAL*8                                      :: ex_fac, tau
  CHARACTER*16                                :: fptoc
  ex_fac = tau/hbar
  DO i=1,n_region
     exp_rt_z(:,i) = evec_region(:,i) * exp ( - eye * eval_region(i) * ex_fac )  
  END DO
  call ecbct(exp_rt,exp_rt_z,evec_region,n_region,n_region,n_region)
  IF(log_main(1)) then
     title='exponential sector propagator at tau = '//fptoc(tau)
     call prntcmn(title,exp_rt,n_region,n_region,n_region,n_region,iout,'e')
  END IF
  END SUBROUTINE umat_reg_h 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck umat_reg_z
!***begin prologue     umat_reg_z
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            regional unitary propagators.
!***
!***references
!***routines called
!***end prologue       umat_reg_z
!
  SUBROUTINE umat_reg_z(eval_region,evec_region_right,     &
                        evec_region_left,exp_rt,exp_rt_z,tau,n_region)
  IMPLICIT NONE
  INTEGER                                     :: n_region, type
  INTEGER                                     :: i
  COMPLEX*16, DIMENSION(n_region)             :: eval_region
  COMPLEX*16, DIMENSION(n_region,n_region)    :: evec_region_left
  COMPLEX*16, DIMENSION(n_region,n_region)    :: evec_region_right
  COMPLEX*16, DIMENSION(n_region,n_region)    :: exp_rt
  COMPLEX*16, DIMENSION(n_region,n_region)    :: exp_rt_z
  REAL*8                                      :: tau, ex_fac
  CHARACTER*16                                :: fptoc
  ex_fac = tau/hbar
  DO i=1,n_region
     exp_rt_z(:,i) = evec_region_right(:,i) * exp ( -eye * eval_region(i) * ex_fac )  
  END DO
!
! Use one of the vectors as scratch.  Its not needed anymore.
!
  call cebhct(exp_rt,exp_rt_z,evec_region_left,n_region,n_region,n_region)
  IF(log_main(1)) then
     title='exponential sector propagator at tau = '//fptoc(tau)
     call prntcmn(title,exp_rt,n_region,n_region,n_region,n_region,iout,'e')
  END IF
  END SUBROUTINE umat_reg_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE regional_module
