!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MODULE dvr_propagators
  USE dvr_global
  USE dvrprop_global
!**begin prologue     dvr_propagators
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords
!**author             schneider, barry (nsf)
!**source
!**purpose            
!**references
!**routines called
!**end prologue       dvr_propagators
!

  INTERFACE diag_reg
     MODULE PROCEDURE diag_reg_d, diag_reg_z
  END INTERFACE diag_reg
  INTERFACE umat_reg
     MODULE PROCEDURE umat_reg_d, umat_reg_z
  END INTERFACE umat_reg
!
  CONTAINS
!
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
!
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
!
!deck diag_reg_d
!***begin prologue     diag_reg_d
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            regional matrix elements of kinetic energy
!                      with zero diagonal is diagonalized.
!***
!***references
!***routines called    dsyev(LAPACK)
!***end prologue       diag_reg_d
!
  SUBROUTINE diag_reg_d(k_region,eval_region,evec_region_right,   &
                        evec_region_left,n_region,reg)
  USE dvr_global
  USE dvrprop_global
  IMPLICIT NONE
  INTEGER                                :: n_region, reg, info
  REAL*8, DIMENSION(n_region)            :: eval_region
  REAL*8, DIMENSION(n_region,n_region)   :: k_region
  REAL*8, DIMENSION(n_region,n_region)   :: evec_region_right
  REAL*8, DIMENSION(n_region,n_region)   :: evec_region_left
  CHARACTER(LEN=4)                       :: itoc
!
!
  evec_region_right = k_region
  call dsyev('v','l',n_region,evec_region_right,n_region,eval_region,   &
              scr_d,5*n_region,info)
  if(log_main(1)) then
     title='sector eigenvalues region = '//itoc(reg)
     call prntfmn(title,eval_region,n_region,1,n_region,1,iout,'e')
     title='sector eigenvectors = '//itoc(reg)
     call prntfmn(title,evec_region_right,n_region,n_region,            &
                                    n_region,n_region,iout,'e')
  end if
  write(iplot(1),*) 'eigenvalues'
  write(iplot(1),*) eval_region
  write(iplot(1),*) 'eigenvectors'
  write(iplot(1),*) evec_region_right
  END SUBROUTINE diag_reg_d  
!
!deck diag_reg_z
!***begin prologue     diag_reg_z
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            regional matrix elements of kinetic energy
!                      with complex potential is diagonalized.
!***
!***references
!***routines called    dsyev(LAPACK)
!***end prologue       diag_reg_z
!
  SUBROUTINE diag_reg_z(k_region,eval_region,evec_region_right, &
                        evec_region_left,n_region,reg)
  USE dvr_global
  USE dvrprop_global
  IMPLICIT NONE
  INTEGER                                  :: n_region, info, reg
  INTEGER                                  :: i, j
  COMPLEX*16, DIMENSION(n_region)          :: eval_region
  COMPLEX*16, DIMENSION(n_region,n_region) :: k_region
  COMPLEX*16, DIMENSION(n_region,n_region) :: evec_region_right
  COMPLEX*16, DIMENSION(n_region,n_region) :: evec_region_left
  COMPLEX*16                               :: cdotc, val_z
!
!
  call zgeev('v','v',n_region,k_region,n_region,eval_region,       &
             evec_region_left,n_region,evec_region_right,n_region, &
             scr_z,10*n_region,scr_z,info)
  IF(info /= 0 ) THEN
     CALL lnkerr('eigenvalue routine error')
  END IF
!
! Renormalize
!
  DO i=1,n_region
     val_z=1.d0/sqrt(cdotc(n_region,evec_region_left(1,i),1,       &
                     evec_region_right(1,i),1))   
     DO j=1,n_region
        evec_region_right(j,i) = val_z * evec_region_right(j,i)
        evec_region_left(j,i) = conjg(val_z) * evec_region_left(j,i)
     END DO
  END DO
  title='sector eigenvalues complex region'
  call prntcmn(title,eval_region,n_region,1,n_region,1,iout,'e')
  if(log_main(1)) then
     title='sector left eigenvectors'
     call prntcmn(title,evec_region_left,n_region,n_region,        &
                  n_region,n_region,iout,'e')
     title='sector right eigenvectors'
     call prntcmn(title,evec_region_right,n_region,n_region,       &
                  n_region,n_region,iout,'e')
  end if
  write(iplot(1),*) 'eigenvalues'
  write(iplot(1),*) eval_region
  write(iplot(1),*) 'left eigenvectors'
  write(iplot(1),*) evec_region_left
  write(iplot(1),*) 'right eigenvectors'
  write(iplot(1),*) evec_region_right
END SUBROUTINE diag_reg_z
END MODULE dvr_propagators
