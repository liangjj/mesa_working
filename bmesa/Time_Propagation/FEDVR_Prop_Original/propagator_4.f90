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
!***references
!***routines called
!***end prologue       propagator_4
!
  SUBROUTINE propagator_4(q)
  USE dvr_global
  USE dvrprop_global
  USE dvr_propagators
  IMPLICIT NONE
  INTEGER                                  :: q, i
  REAL*8, DIMENSION(4)                     :: tau
  REAL*8                                   :: p
!
  p = 1.d0/( 4.d0 - 4.d0**(1.d0/3.d0) )
!
!
  DO i=1,num_reg(q),2
     ALLOCATE(mat_reg_d(i,q)%cosine_t_mat(nfun_reg(i,q),nfun_reg(i,q),3), &
              mat_reg_d(i,q)%sine_t_mat(nfun_reg(i,q),nfun_reg(i,q),3))
  END DO
  DO i=2,num_reg(q),2
     ALLOCATE(mat_reg_d(i,q)%cosine_t_mat(nfun_reg(i,q),nfun_reg(i,q),2), &
              mat_reg_d(i,q)%sine_t_mat(nfun_reg(i,q),nfun_reg(i,q),2))
  END DO
  ALLOCATE(si_d(maxmem_reg,maxmem_reg),ci_d(maxmem_reg,maxmem_reg))
!
! Construct the propagators
!
  tau(1) = p * deltat *.5d0
  tau(2) = p * deltat
  tau(3) = ( 1.d0 - 4.d0 * p ) * deltat *.5d0
  tau(4) = ( 1.d0 - 4.d0 * p ) * deltat
  DO i=1,num_reg(q),2
     call umat_reg(mat_reg_d(i,q)%eigval_mat_d,          &
                   mat_reg_d(i,q)%eigvec_mat_d,          &
                   mat_reg_d(i,q)%eigvec_mat_d,          &
                   mat_reg_d(i,q)%cosine_t_mat(1,1,1),   &
                   mat_reg_d(i,q)%sine_t_mat(1,1,1),     &
                   si_d,ci_d,tau(1),nfun_reg(i,q),i)
     call umat_reg(mat_reg_d(i,q)%eigval_mat_d,          &
                   mat_reg_d(i,q)%eigvec_mat_d,          &
                   mat_reg_d(i,q)%eigvec_mat_d,          &
                   mat_reg_d(i,q)%cosine_t_mat(1,1,2),   &
                   mat_reg_d(i,q)%sine_t_mat(1,1,2),     &
                   si_d,ci_d,tau(2),nfun_reg(i,q),i)
     call umat_reg(mat_reg_d(i,q)%eigval_mat_d,          &
                   mat_reg_d(i,q)%eigvec_mat_d,          &
                   mat_reg_d(i,q)%eigvec_mat_d,          &
                   mat_reg_d(i,q)%cosine_t_mat(1,1,3),   &
                   mat_reg_d(i,q)%sine_t_mat(1,1,3),     &
                   si_d,ci_d,tau(4),nfun_reg(i,q),i)
  END DO
  DO i=2,num_reg(q),2
     call umat_reg(mat_reg_d(i,q)%eigval_mat_d,          &
                   mat_reg_d(i,q)%eigvec_mat_d,          &
                   mat_reg_d(i,q)%eigvec_mat_d,          &
                   mat_reg_d(i,q)%cosine_t_mat(1,1,1),   &
                   mat_reg_d(i,q)%sine_t_mat(1,1,1),     &
                   si_d,ci_d,tau(2),nfun_reg(i,q),i)
     call umat_reg(mat_reg_d(i,q)%eigval_mat_d,          &
                   mat_reg_d(i,q)%eigvec_mat_d,          &
                   mat_reg_d(i,q)%eigvec_mat_d,          &
                   mat_reg_d(i,q)%cosine_t_mat(1,1,2),   &
                   mat_reg_d(i,q)%sine_t_mat,            &
                   si_d,ci_d,tau(4),nfun_reg(i,q),i)
  END DO
  DEALLOCATE(si_d,ci_d)
END SUBROUTINE propagator_4
