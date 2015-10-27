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
!deck propagator_2.f
!***begin prologue     propagator_2
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            second order propagators
!***
!***references
!***routines called
!***end prologue       propagator_2
!
  SUBROUTINE propagator_2(q)
  USE dvr_global
  USE dvrprop_global
  USE regional_umat
  IMPLICIT NONE
  INTEGER                                  :: q, i, even_odd
  REAL*8                                   :: tau
!
!
! Construct the "outer" propagator at delta_t/2
!
  even_odd = num_reg(q) - 2 * ( num_reg(q) / 2 )
  IF (even_odd == 0 ) THEN
      tau = deltat*.5d0
      DO i=1,num_reg(q),2
         call umat_reg(mat_reg_d(i,q)%eigval_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%cosine_t_mat,      &
                       mat_reg_d(i,q)%sine_t_mat,        &
                       si_d,ci_d,tau,nfun_reg(i,q),1)
      END DO
!
! Construct the "inner" propagator at delta_t
!
      tau = deltat
      DO i=2,num_reg(q)-1,2
         call umat_reg(mat_reg_d(i,q)%eigval_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%cosine_t_mat,      &
                       mat_reg_d(i,q)%sine_t_mat,        &
                       si_d,ci_d,tau,nfun_reg(i,q),1)
      END DO
      i = num_reg(q)
      IF(.not.absorb) THEN
         call umat_reg(mat_reg_d(i,q)%eigval_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%cosine_t_mat,      &
                       mat_reg_d(i,q)%sine_t_mat,        &
                       si_d,ci_d,tau,nfun_reg(i,q),1)
      ELSE
         call umat_reg(mat_reg_z(q)%eigval_mat_z,        &
                       mat_reg_z(q)%eigvec_mat_z_r,      &
                       mat_reg_z(q)%eigvec_mat_z_l,      &
                       mat_reg_d(i,q)%cosine_t_mat,      &
                       mat_reg_d(i,q)%sine_t_mat,        &
                       si_z,ci_z,tau,nfun_reg(i,q),1)
      END IF
  ELSE
      tau = deltat*.5d0
      DO i=1,num_reg(q)-1,2
         call umat_reg(mat_reg_d(i,q)%eigval_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%eigvec_mat_d,      &
                       mat_reg_d(i,q)%cosine_t_mat,      &
                       mat_reg_d(i,q)%sine_t_mat,        &
                       si_d,ci_d,tau,nfun_reg(i,q),1)
      END DO
      i = num_reg(q)
      IF(.not.absorb) THEN
          call umat_reg(mat_reg_d(i,q)%eigval_mat_d,     &
                        mat_reg_d(i,q)%eigvec_mat_d,     &
                        mat_reg_d(i,q)%eigvec_mat_d,     &
                        mat_reg_d(i,q)%cosine_t_mat,     &
                        mat_reg_d(i,q)%sine_t_mat,       &
                        si_d,ci_d,tau,nfun_reg(i,q),1)
      ELSE
          call umat_reg(mat_reg_z(q)%eigval_mat_z,       &
                        mat_reg_z(q)%eigvec_mat_z_r,     &
                        mat_reg_z(q)%eigvec_mat_z_l,     &
                        mat_reg_d(i,q)%cosine_t_mat,     &
                        mat_reg_d(i,q)%sine_t_mat,       &
                        si_z,ci_z,tau,nfun_reg(i,q),1)
      END IF
!
! Construct the "inner" propagator at delta_t
!
      tau = deltat
      DO i=2,num_reg(q),2
         call umat_reg(mat_reg_d(i,q)%eigval_mat_d,      &
                         mat_reg_d(i,q)%eigvec_mat_d,    &
                         mat_reg_d(i,q)%eigvec_mat_d,    &
                         mat_reg_d(i,q)%cosine_t_mat,    &
                         mat_reg_d(i,q)%sine_t_mat,      &
                         si_d,ci_d,tau,nfun_reg(i,q),1)
      END DO
  END IF
END SUBROUTINE propagator_2
