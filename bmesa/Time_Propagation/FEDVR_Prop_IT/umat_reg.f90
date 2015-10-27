!***begin prologue     umat_reg
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            regional unitary propagators.
!***
!***references
!***routines called
!***end prologue       umat_reg
!
  SUBROUTINE umat_reg(eval_region,evec_region,exp_t,exp_t_d,tau,n_region)
  USE dvrprop_global_it
  IMPLICIT NONE
  INTEGER                                     :: n_region
  INTEGER                                     :: i
  REAL*8, DIMENSION(n_region)                 :: eval_region
  REAL*8, DIMENSION(n_region,n_region)        :: evec_region
  REAL*8, DIMENSION(n_region,n_region)        :: exp_t
  REAL*8, DIMENSION(n_region,n_region)        :: exp_t_d
  REAL*8                                      :: ex_fac, tau
  CHARACTER*16                                :: fptoc
  ex_fac = tau/hbar
  do i=1,n_region
     exp_t_d(:,i) = evec_region(:,i) * exp(- eval_region(i)*ex_fac)  
  end do
  call ebct(exp_t,exp_t_d,evec_region,n_region,n_region,n_region)
  if(log_main(1)) then
     title='exponential sector propagator at tau = '//fptoc(tau)
     call prntfmn(title,exp_t,n_region,n_region,                 &
                  n_region,n_region,iout,'e')
  end if
  END SUBROUTINE umat_reg 
