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
  IMPLICIT NONE
  INTEGER                                  :: q, i, even_odd
  REAL*8                                   :: tau
  REAL*8, DIMENSION(:,:), ALLOCATABLE      :: si, ci 
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE  :: si_c, ci_c 
!
!
  DO i=1,num_reg(q)
     ALLOCATE(                                                   &
     mat_reg_d(i,q)%cosine_t_mat(nfun_reg(i,q),nfun_reg(i,q),1), &
     mat_reg_d(i,q)%sine_t_mat(nfun_reg(i,q),nfun_reg(i,q),1))
  END DO
  ALLOCATE(si_d(maxmem_reg,maxmem_reg),                     &
           ci_d(maxmem_reg,maxmem_reg))
  IF(absorb) THEN
     ALLOCATE(si_z(maxmem_reg,maxmem_reg),                  &
              ci_z(maxmem_reg,maxmem_reg))
  END IF
!
!
! Construct the "outer" propagator at delta_t/2
!
  tau = deltat*.5d0
  DO i=1,num_reg(q),2
     call umat_reg_d(mat_reg_d(i,q)%eigval_mat_d,           &
                     mat_reg_d(i,q)%eigvec_mat_d,           &
                     mat_reg_d(i,q)%eigvec_mat_d,           &
                     mat_reg_d(i,q)%cosine_t_mat,         &
                     mat_reg_d(i,q)%sine_t_mat,           &
                     si_d,ci_d,tau,nfun_reg(i,q))
  END DO
!
! Construct the "inner" propagator at delta_t
!
  tau = deltat
  DO i=2,num_reg(q),2
     call umat_reg_d(mat_reg_d(i,q)%eigval_mat_d,           &
                   mat_reg_d(i,q)%eigvec_mat_d,             &
                   mat_reg_d(i,q)%eigvec_mat_d,             &
                   mat_reg_d(i,q)%cosine_t_mat,           &
                   mat_reg_d(i,q)%sine_t_mat,             &
                   si_d,ci_d,tau,nfun_reg(i,q))
  END DO
  tau = deltat*.5d0
  DO i=1,num_reg(q),2
     call umat_reg_d(mat_reg_d(i,q)%eigval_mat_d,           &
                   mat_reg_d(i,q)%eigvec_mat_d,             &
                   mat_reg_d(i,q)%eigvec_mat_d,             &
                   mat_reg_d(i,q)%cosine_t_mat,           &
                   mat_reg_d(i,q)%sine_t_mat,               &
                   si_d,ci_d,tau,nfun_reg(i,q))
  END DO
  DEALLOCATE(si_d,ci_d)
  IF(absorb) THEN
     DEALLOCATE(si_z,ci_z)
  END IF
END SUBROUTINE propagator_2
