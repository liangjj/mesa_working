!***********************************************************************
! scattering
!**begin prologue     scattering
!**date written       082805   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Time_Propagation
!**purpose            Contains all of the subroutines to compute
!***                  scattering phase shifts.
!***                  Explicit interfaces are used to allow
!***                  a transparent use of generic subroutines which work
!***                  for both real and complex vectors.  This feature
!***                  permits a single code to be used for both real and
!***                  imaginary time propagation.
!***description       See subroutines
!**references
!**modules needed     See USE statements below
!**end prologue       scattering_module
!***********************************************************************
!deck scattering.f
!***begin prologue     scattering
!***date written       040706   (yymmdd)
!***revision date               (yymmdd)
!***keywords           space, propagation, dvr, finite difference
!***author             schneider, b. i.(nsf)
!***source
!***purpose            the global FEDVR matrices are split into their
!***                   regional parts for ease of use in the propagation.
!***                   the regional kinetic energy matrix may be modified
!***                   to include any diagonal contribution from the one
!***                   body potential or any other diagnal part.  once this
!***                   is done, the kinetic energy is diagonalized.  the 
!***                   eigenvectors and eigenvalues are used to form the
!***                   element propagators later in the code.
!***references
!***routines called    iosys, util and mdutil
!***end prologue       space
  Subroutine scattering (eig_val,surface_vec,k,l_val,potential_type,n)
  USE dvr_global
  USE dvr_shared
  USE dvrprop_global
  IMPLICIT NONE
  REAL*8, DIMENSION(:)                     :: eig_val
  REAL*8, DIMENSION(:)                     :: surface_vec
  REAL*8                                   :: k
  REAL*8                                   :: energy
  REAL*8                                   :: r_matrix
  INTEGER                                  :: l_val
  INTEGER                                  :: i, j, start, row_dim
  CHARACTER(LEN=*)                         :: potential_type
!
!
  energy = k*k*.5d0
  r_matrix = 0.0d0
  DO i=1,n_eig
     r_matrix = r_matrix + surface_vec(i) * surface_vec(i)   &
                                          /                  &
                             ( energy - eig_val(i) )
  END DO
  IF ( type_potential /= 'coulomb' ) then

  ELSE 

1 FORMAT(/15x,'Regional Information for DVR Basis,'                &
              ' Variable = ',i2)
2 FORMAT(/15x,'Regional Information for FD Basis,'                 &
                 ' Variable = ',i2)
3 FORMAT(/,17x,'Region',5x,'Number of Functions',5x,'Type')
4 FORMAT(16x,i4,14x,i5,9x,a8)
END Subroutine scattering
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!deck regional_diag_d
!***begin prologue     regional_diag_d
!***date written       040706   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            Regional matrix elements of kinetic energy
!                      with zero diagonal are diagonalized.
!***
!***references
!***routines called    dsyev(LAPACK)
!***end prologue       regional_diag_d
!
  SUBROUTINE regional_diag_d(k_region,eval_region,evec_region,         &
                           n_region,reg)
  IMPLICIT NONE
  INTEGER                                :: n_region, reg
  INTEGER                                :: info
  REAL*8, DIMENSION(n_region,n_region)   :: k_region
  REAL*8, DIMENSION(n_region)            :: eval_region
  REAL*8, DIMENSION(n_region,n_region)   :: evec_region
  CHARACTER(LEN=4)                       :: itoc
!
!
  evec_region = k_region
  call dsyev('v','l',n_region,evec_region,n_region,eval_region,   &
              scr_d,5*n_region,info)
  if(log_main(1)) then
     title='sector eigenvalues region = '//itoc(reg)
     call prntfmn(title,eval_region,n_region,1,n_region,1,iout,'e')
     title='sector eigenvectors = '//itoc(reg)
     call prntfmn(title,evec_region,n_region,n_region,            &
                                    n_region,n_region,iout,'e')
  end if
  END SUBROUTINE regional_diag_d  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck regional_diag_z
!***begin prologue     regional_diag_z
!***date written       040706   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            regional matrix elements of kinetic energy
!                      with complex potential is diagonalized.
!***
!***references
!***routines called    dsyev(LAPACK)
!***end prologue       regional_diag_z
!
  SUBROUTINE regional_diag_z(k_region,eval_region,evec_region_right,    &
                             evec_region_left,n_region,reg)
  IMPLICIT NONE
  INTEGER                                  :: n_region, reg
  INTEGER                                  :: info
  INTEGER                                  :: i, j
  COMPLEX*16, DIMENSION(n_region,n_region) :: k_region
  COMPLEX*16, DIMENSION(n_region)          :: eval_region
  COMPLEX*16, DIMENSION(n_region,n_region) :: evec_region_right
  COMPLEX*16, DIMENSION(n_region,n_region) :: evec_region_left
  COMPLEX*16                               :: val_z, cdotc
!
!
  call zgeev('v','v',n_region,k_region,n_region,eval_region,       &
             evec_region_left,n_region,evec_region_right,n_region, &
             scr_z,10*n_region,scr_d,info)
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
END SUBROUTINE regional_diag_z
!
!***********************************************************************
!***********************************************************************

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE regional_module
