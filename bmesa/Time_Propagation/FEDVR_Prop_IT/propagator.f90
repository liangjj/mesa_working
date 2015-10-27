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
  SUBROUTINE propagator(q)
  USE dvr_global
  USE dvrprop_global_it
  IMPLICIT NONE
  INTEGER                                  :: q, i, even_odd
  INTEGER                                  :: trips, ntrip
  REAL*8                                   :: tau
!
!
! Special case of only one region
!
  write(iout,1) prop_order, p_fac
!
  IF ( num_reg(q) == 1) THEN
       write(iout,2)
       IF(prop_order==2) THEN
          tau = deltat
          call umat_reg(mat_reg_d(1,q)%eigval_mat_d,             &
                        mat_reg_d(1,q)%eigvec_mat_d,             &
                        mat_reg_d(1,q)%exp_t_mat(:,:,1),         &
                        exp_tmp_d,tau,nfun_reg(1,q))
       ELSE IF(prop_order==4) THEN
          p_loc = p_fac
          trips = 2
          Write(iout,3) trips
          DO ntrip=1,trips
             write(iout,4) ntrip, p_loc
              tau = p_loc * deltat
              call umat_reg(mat_reg_d(1,q)%eigval_mat_d,          &
                            mat_reg_d(1,q)%eigvec_mat_d,          &
                            mat_reg_d(1,q)%exp_t_mat(:,:,ntrip),  &
                            exp_tmp_d,tau,nfun_reg(1,q))
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
          tau = p_loc * deltat *.5d0
          DO i=1,num_reg(q),2
             call umat_reg(mat_reg_d(i,q)%eigval_mat_d,             &
                           mat_reg_d(i,q)%eigvec_mat_d,             &
                           mat_reg_d(i,q)%exp_t_mat(:,:,ntrip),     &
                           exp_tmp_d,tau,nfun_reg(i,q))
          END DO
!
!         Construct the "inner" propagator at p_loc * delta_t
!
          tau = p_loc * deltat
          DO i=2,num_reg(q),2
             call umat_reg(mat_reg_d(i,q)%eigval_mat_d,             &
                           mat_reg_d(i,q)%eigvec_mat_d,             &
                           mat_reg_d(i,q)%exp_t_mat(:,:,ntrip),     &
                           exp_tmp_d,tau,nfun_reg(i,q))
          END DO
!
!         Modify p_loc for fourth order propagator so ist correct second time 
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
END SUBROUTINE propagator
