!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             MODULE regional_umat
                             USE dvrprop_global
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
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
END MODULE regional_umat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
