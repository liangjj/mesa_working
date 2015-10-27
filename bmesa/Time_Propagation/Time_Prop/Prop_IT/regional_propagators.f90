!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          MODULE regional_propagators
                               USE dvr_global
                               USE dvrprop_global
                               USE regional_umat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck reg_prop_d.f
!***begin prologue     reg_prop_d
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
  SUBROUTINE reg_prop_d(q,t_in)
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
END SUBROUTINE reg_prop_d
!***********************************************************************
!***********************************************************************
!deck reg_prop_h.f
!***begin prologue     reg_prop_h
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
  SUBROUTINE reg_prop_h(q,t_in)
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
END SUBROUTINE reg_prop_h
!
!***********************************************************************
!***********************************************************************
  SUBROUTINE reg_prop_z(q,t_in)
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
END SUBROUTINE reg_prop_z
!
!***********************************************************************
!***********************************************************************
                   END MODULE regional_propagators
