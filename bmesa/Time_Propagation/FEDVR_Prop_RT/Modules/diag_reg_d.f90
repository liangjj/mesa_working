!deck diag_reg_d
!***begin prologue     diag_reg_d
!***date written       040706   (yymmdd)
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
  SUBROUTINE diag_reg_d(k_region,eval_region,evec_region_right,         &
                        evec_region_left,n_region,reg)
  USE dvr_global
  USE dvrprop_global_rt
  IMPLICIT NONE
  INTEGER                                :: n_region, reg
  INTEGER                                :: info
  REAL*8, DIMENSION(n_region,n_region)   :: k_region
  REAL*8, DIMENSION(n_region)            :: eval_region
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
