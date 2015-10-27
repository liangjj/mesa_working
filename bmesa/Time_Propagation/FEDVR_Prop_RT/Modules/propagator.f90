! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{Unitary Time Propagators}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck propagator.f
!***begin prologue     propagator
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            propagators
!***
!***description        The second order decomposition of the time exponential
!***                   is performed for a Hamiltonian H = H_a + H_b + V as
!***                   U_2(tau) = exp(-iV * tau/2)   *
!***                                    u_2(tau)     *
!***                              exp(-iV * tau/2)
!***
!***                   u_2(tau) = exp(-iH_a * tau/2) * 
!***                              exp(-iH_b * tau)   * 
!***                              exp(_iH_a * tau/2) *
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
  SUBROUTINE propagator(q)
  USE dvr_global
  USE dvrprop_global_rt
  USE regional_umat
  IMPLICIT NONE
  INTEGER                                  :: q, i, even_odd
  INTEGER                                  :: trips, ntrip
  REAL*8                                   :: tau
!
!
! Special case of only one region
!
  IF ( num_reg(q) == 1) THEN
       tau = deltat
       call umat_reg(mat_reg_d(1,q)%eigval_mat_d,             &
                     mat_reg_d(1,q)%eigvec_mat_d,             &
                     mat_reg_d(1,q)%eigvec_mat_d,             &
                     mat_reg_d(1,q)%cosine_t_mat(:,:,1),      &
                     mat_reg_d(1,q)%sine_t_mat(:,:,1),        &
                     si_d,ci_d,tau,nfun_reg(1,q))
  ELSE
       p_loc = p_fac
       trips = 1
       IF ( prop_order == 4 ) THEN
            trips = 2
       END IF
!
       DO ntrip=1,trips
!
!         Construct the "outer" propagator at p_fac * delta_t/2
!
          even_odd = num_reg(q) - 2 * ( num_reg(q) / 2 )
          IF (even_odd == 0 ) THEN
              tau = p_loc * deltat *.5d0
              DO i=1,num_reg(q),2
                 call umat_reg(mat_reg_d(i,q)%eigval_mat_d,             &
                               mat_reg_d(i,q)%eigvec_mat_d,             &
                               mat_reg_d(i,q)%eigvec_mat_d,             &
                               mat_reg_d(i,q)%cosine_t_mat(:,:,ntrip),  &
                               mat_reg_d(i,q)%sine_t_mat(:,:,ntrip),    &
                               si_d,ci_d,tau,nfun_reg(i,q))
              END DO
!
!         Construct the "inner" propagator at p_loc * delta_t
!
              tau = p_loc * deltat
              DO i=2,num_reg(q)-1,2
                 call umat_reg(mat_reg_d(i,q)%eigval_mat_d,             &
                               mat_reg_d(i,q)%eigvec_mat_d,             &
                               mat_reg_d(i,q)%eigvec_mat_d,             &
                               mat_reg_d(i,q)%cosine_t_mat(:,:,ntrip),  &
                               mat_reg_d(i,q)%sine_t_mat(:,:,ntrip),    &
                               si_d,ci_d,tau,nfun_reg(i,q))
              END DO
              i = num_reg(q)
              IF(.not.absorb) THEN
                 call umat_reg(mat_reg_d(i,q)%eigval_mat_d,             &
                               mat_reg_d(i,q)%eigvec_mat_d,             &
                               mat_reg_d(i,q)%eigvec_mat_d,             &
                               mat_reg_d(i,q)%cosine_t_mat(:,:,ntrip),  &
                               mat_reg_d(i,q)%sine_t_mat(:,:,ntrip),    &
                               si_d,ci_d,tau,nfun_reg(i,q))
              ELSE
                 call umat_reg(mat_reg_z(q)%eigval_mat_z,               &
                               mat_reg_z(q)%eigvec_mat_z_r,             &
                               mat_reg_z(q)%eigvec_mat_z_l,             &
                               mat_reg_d(i,q)%cosine_t_mat(:,:,ntrip),  &
                               mat_reg_d(i,q)%sine_t_mat(:,:,ntrip),    &
                              si_z,ci_z,tau,nfun_reg(i,q))
              END IF
          ELSE
              tau = p_loc * deltat *.5d0
              DO i=1,num_reg(q)-1,2
                 call umat_reg(mat_reg_d(i,q)%eigval_mat_d,             &
                               mat_reg_d(i,q)%eigvec_mat_d,             &
                               mat_reg_d(i,q)%eigvec_mat_d,             &
                               mat_reg_d(i,q)%cosine_t_mat(:,:,ntrip),  &
                               mat_reg_d(i,q)%sine_t_mat(:,:,ntrip),    &
                               si_d,ci_d,tau,nfun_reg(i,q))
              END DO
              i = num_reg(q)
              IF(.not.absorb) THEN
                 call umat_reg(mat_reg_d(i,q)%eigval_mat_d,             &
                               mat_reg_d(i,q)%eigvec_mat_d,             &
                               mat_reg_d(i,q)%eigvec_mat_d,             &
                               mat_reg_d(i,q)%cosine_t_mat(:,:,ntrip),  &
                               mat_reg_d(i,q)%sine_t_mat(:,:,ntrip),    &
                               si_d,ci_d,tau,nfun_reg(i,q))
              ELSE
                 call umat_reg(mat_reg_z(q)%eigval_mat_z,               &
                               mat_reg_z(q)%eigvec_mat_z_r,             &
                               mat_reg_z(q)%eigvec_mat_z_l,             &
                               mat_reg_d(i,q)%cosine_t_mat(:,:,ntrip),  &
                               mat_reg_d(i,q)%sine_t_mat(:,:,ntrip),    &
                              si_z,ci_z,tau,nfun_reg(i,q))
              END IF
!
!         Construct the "inner" propagator at delta_t
!
              tau = p_loc * deltat
              DO i=2,num_reg(q),2
                 call umat_reg(mat_reg_d(i,q)%eigval_mat_d,             &
                               mat_reg_d(i,q)%eigvec_mat_d,             &
                               mat_reg_d(i,q)%eigvec_mat_d,             &
                               mat_reg_d(i,q)%cosine_t_mat(:,:,ntrip),  &
                               mat_reg_d(i,q)%sine_t_mat(:,:,ntrip),    &
                               si_d,ci_d,tau,nfun_reg(i,q))
              END DO
          END IF
!
!         Modify p_loc for fourth order propagator so ist correct second time 
!         through last loop.
!
          p_loc = ( 1.d0 - 4.d0 * p_loc )
       END DO
  END IF
END SUBROUTINE propagator
