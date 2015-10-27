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
  write(iplot(1),*) 'eigenvalues'
  write(iplot(1),*) eval_region
  write(iplot(1),*) 'left eigenvectors'
  write(iplot(1),*) evec_region_left
  write(iplot(1),*) 'right eigenvectors'
  write(iplot(1),*) evec_region_right
END SUBROUTINE diag_reg_z