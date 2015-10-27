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
!deck propagator_4.f
!***begin prologue     propagator_4
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            fourth order propagators
!***
!***description        The fourth order decomposition of the time exponential
!***                   is performed for a Hamiltonian H = H_a + H_b as
!***                   U_4(tau) = U_2(p*tau) * U_2(p*tau) * U_2((1-4*p)*tau) * U_2(p*tau) * U_2(p*tau) 
!***                                       Where
!***                   U_2(tau) =  exp(-iV * tau/2)   * 
!***                               exp(-iH_a * tau/2) * 
!***                               exp(-iH_b * tau)   * 
!***                               exp(-iH_a * tau/2) *
!***                               exp(-iV * uau/2)
!***                            p = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
!***
!***                   The expression for U_4 requires the calculation of four 
!***                   independent propagators.
!***references
!***routines called
!***end prologue       propagator_4
!
  SUBROUTINE propagator_4(q)
  USE dvr_global
  USE dvrprop_global_rt
  USE regional_umat
  IMPLICIT NONE
  INTEGER                                  :: q, i
  REAL*8, DIMENSION(4)                     :: tau
  REAL*8                                   :: p
!
  IF(num_reg(q) == 1 ) THEN
     tau(1) = deltat
     call umat_reg(mat_reg_d(1,q)%eigval_mat_d,            &
                   mat_reg_d(1,q)%eigvec_mat_d,            &
                   mat_reg_d(1,q)%eigvec_mat_d,            &
                   mat_reg_d(1,q)%cosine_t_mat(:,:,1),     &
                   mat_reg_d(1,q)%sine_t_mat(:,:,1),       &
                   si_d,ci_d,tau(1),nfun_reg(1,q))
  ELSE
!
     p = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
!
!
!
! Construct the propagators
!
     tau(1) = p * deltat *.5d0
     tau(2) = p * deltat
!
!
     tau(3) = ( 1.d0 - 4.d0 * p ) * deltat *.5d0
     tau(4) = ( 1.d0 - 4.d0 * p ) * deltat
     DO i=1,num_reg(q),2
!
!       Odd propagator at p * deltat * .5d0
!
        call umat_reg(mat_reg_d(i,q)%eigval_mat_d,            &
                      mat_reg_d(i,q)%eigvec_mat_d,            &
                      mat_reg_d(i,q)%eigvec_mat_d,            &
                      mat_reg_d(i,q)%cosine_t_mat(:,:,1),     &
                      mat_reg_d(i,q)%sine_t_mat(:,:,1),       &
                      si_d,ci_d,tau(1),nfun_reg(i,q))
!
!       Odd propagator at ( 1.d0 - 4 * p ) * deltat * .5d0
!
        call umat_reg(mat_reg_d(i,q)%eigval_mat_d,            &
                      mat_reg_d(i,q)%eigvec_mat_d,            &
                      mat_reg_d(i,q)%eigvec_mat_d,            &
                      mat_reg_d(i,q)%cosine_t_mat(:,:,2),     &
                      mat_reg_d(i,q)%sine_t_mat(:,:,2),       &
                      si_d,ci_d,tau(3),nfun_reg(i,q))
     END DO
     DO i=2,num_reg(q),2
!
!       Even propagator at  p * deltat
!
        call umat_reg(mat_reg_d(i,q)%eigval_mat_d,            &
                      mat_reg_d(i,q)%eigvec_mat_d,            &
                      mat_reg_d(i,q)%eigvec_mat_d,            &
                      mat_reg_d(i,q)%cosine_t_mat(:,:,1),     &
                      mat_reg_d(i,q)%sine_t_mat(:,:,1),       &
                      si_d,ci_d,tau(2),nfun_reg(i,q))
!
!       Even propagator at  ( 1.d0 - 4.d0 * p ) * deltat
!
        call umat_reg(mat_reg_d(i,q)%eigval_mat_d,            &
                      mat_reg_d(i,q)%eigvec_mat_d,            &
                      mat_reg_d(i,q)%eigvec_mat_d,            &
                      mat_reg_d(i,q)%cosine_t_mat(:,:,2),     &
                      mat_reg_d(i,q)%sine_t_mat(:,:,2),       &
                      si_d,ci_d,tau(4),nfun_reg(i,q))
     END DO
  END IF
END SUBROUTINE propagator_4
