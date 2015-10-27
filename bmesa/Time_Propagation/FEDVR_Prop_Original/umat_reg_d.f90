!deck umat_reg_d
!***begin prologue     umat_reg_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            regional unitary propagators.
!***
!***references
!***routines called
!***end prologue       umat_reg_d
!
  SUBROUTINE umat_reg_d(eval_region,evec_region_right,evec_region_left,   &
                        cos_t,sin_t,si,ci,tau,n_region,reg)
  USE dvr_global
  USE dvrprop_global
  IMPLICIT NONE
  INTEGER                                :: n_region, reg, i, j
  REAL*8, DIMENSION(n_region)            :: eval_region
  REAL*8, DIMENSION(n_region,n_region)   :: evec_region_right
  REAL*8, DIMENSION(n_region,n_region)   :: evec_region_left
  REAL*8, DIMENSION(n_region,n_region)   :: cos_t, sin_t
  REAL*8, DIMENSION(n_region,n_region)   :: si, ci
  REAL*8                                 :: tau, ex_fac
  CHARACTER*16                           :: fptoc
!
!
  ex_fac = tau/hbar
  do i=1,n_region
     ci(:,i) =    evec_region_right(:,i) * cos(eval_region(i)*ex_fac)  
     si(:,i) =  - evec_region_right(:,i) * sin(eval_region(i)*ex_fac)  
  end do
  call ebct(cos_t,ci,evec_region_right,n_region,n_region,n_region)
  call ebct(sin_t,si,evec_region_right,n_region,n_region,n_region)
  if(log_main(1)) then
     title='cosine sector propagator at tau = '//fptoc(tau)
     call prntfmn(title,cos_t,n_region,n_region,                 &
                  n_region,n_region,iout,'e')
     title='sine sector propagator at tau = '//fptoc(tau)
     call prntfmn(title,sin_t,n_region,n_region,                 &
                  n_region,n_region,iout,'e')
  end if
  END SUBROUTINE umat_reg_d 
