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
  SUBROUTINE umat_reg_z(eval_region,evec_region_right,         &
                        evec_region_left,cos_t,sin_t,          &
                        si,ci,tau,n_region,reg)
  USE dvr_global
  USE dvrprop_global
  IMPLICIT NONE
  INTEGER                                   :: n_region, reg, i, j
  COMPLEX*16, DIMENSION(n_region)           :: eval_region
  COMPLEX*16, DIMENSION(n_region,n_region)  :: evec_region_left
  COMPLEX*16, DIMENSION(n_region,n_region)  :: evec_region_right
  REAL*8,     DIMENSION(n_region,n_region)  :: cos_t, sin_t
  COMPLEX*16, DIMENSION(n_region,n_region)  :: si, ci
  REAL*8                                    :: tau, ex_fac
  REAL*8                                    :: re_eig, im_eig
  CHARACTER*16                              :: fptoc
!
!
  call cehbtc(ci,evec_region_left,evec_region_right,         &
              n_region,n_region,n_region)
  call cebhct(ci,evec_region_right,evec_region_left,         &
              n_region,n_region,n_region)
  ex_fac = tau/hbar
  do i=1,n_region
     re_eig = real(eval_region(i))
     im_eig = exp ( imag(eval_region(i)) )
     ci(:,i) =    im_eig * evec_region_right(:,i)              &
                         *                                     &
                  cos(re_eig*ex_fac)  
     si(:,i) =  - im_eig * evec_region_right(:,i)              &
                         *                                     &
                  sin(re_eig*ex_fac)  
  end do
!
! Use one of the vectors as scratch.  Its not needed anymore.
!
  call cebhct(evec_region_right,ci,evec_region_left,            &
              n_region,n_region,n_region)
  cos_t(:,:) = real(evec_region_right(:,:))
  call cebhct(evec_region_right,si,evec_region_left,            &
              n_region,n_region,n_region)
  sin_t(:,:) = imag(evec_region_right(:,:))
  if(log_main(1)) then
     title='cosine sector propagator at tau = '//fptoc(tau)
     call prntfmn(title,cos_t,n_region,n_region,                 &
                  n_region,n_region,iout,'e')
     title='sine sector propagator at tau = '//fptoc(tau)
     call prntfmn(title,sin_t,n_region,n_region,                 &
                  n_region,n_region,iout,'e')
  end if
  END SUBROUTINE umat_reg_z
