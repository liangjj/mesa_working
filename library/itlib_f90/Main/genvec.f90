!deck genvec.f
!***begin prologue     genvec
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           davidson, trial, vector
!***author             schneider, barry (nsf)
!***source
!***purpose            generate a new trial vector based on some zeroth
!***                   order model.  this is equivalent to a preconditioning
!***                   of the matrix.
!***references

!***routines called
!***end prologue       genvec
  SUBROUTINE genvec(vin)
  USE io
  USE dvr_shared
  USE dvd_global
  USE dvd_prnt
  USE precond_global
  IMPLICIT NONE
  REAL*8                                   :: test
  CHARACTER (LEN=4)                        :: itoc
  REAL*8, DIMENSION(n_dvd,addvec)          :: vin
  IF(log_dvd(10)) THEN
     title='input residuals to genvec'
     CALL prntrm(title,resid,n_dvd,addvec,n_dvd,addvec,iout)
  END IF 
  IF(precon == 'diagonal') THEN
     DO  i = 1, addvec
         DO  j=1,n_dvd
            test=eig(i) - diag(j)
            IF(ABS(test) >= nrzero) THEN
               vin(j,i) = resid(j,i)/test
            ELSE
               vin(j,i) = one
            END IF
         END DO
     END DO
     title='vectors'
     CALL prntrm(title,vin,n_dvd,addvec,n_dvd,addvec,iout)
  ELSE IF(precon == 'separable') THEN
!
!    The idea is to take the residual from the DVR basis to
!    something better based on a separable model.  So, the
!    residual get transformed, the result is divided by the
!    diagonal of the separable model and then that is transformed
!    back to the DVR basis.
!
     title='one-body eigenvalues'
     CALL prntrm(title,grid(1)%eigv,nphy(1),1,nphy(1),1,iout)
     title='eigenvalues'
     CALL prntrm(title,eig,addvec,1,addvec,1,iout)
     title='one-body transformation matrix'
     CALL prntrm(title,grid(1)%eigvec,nphy(1),nphy(1),nphy(1), &
                                              nphy(1),iout)
     title='residuals'
     CALL prntrm(title,resid,nphy(1),addvec,nphy(1),addvec,iout)
     IF(dim == 1) THEN
        CALL h12h0(resid,work)
        title='h12h0 work'
        CALL prntrm(title,work,nphy(1),addvec,nphy(1),addvec,iout)
        CALL h1e(work)
        title='h1e work'
        CALL prntrm(title,work,nphy(1),addvec,nphy(1),addvec,iout)
        CALL h12dvr(vin,work)
        title='vectors'
        CALL prntrm(title,vin,nphy(1),addvec,nphy(1),addvec,iout)
     ELSE IF(dim == 2) THEN
        CALL h22h0(resid,resid,work)
        CALL h2e(resid)
        CALL h22dvr(vin,resid,work)
     ELSE IF(dim == 3) THEN
        CALL h32h0(resid,resid,work,svec)
        CALL h3e(resid)
        CALL h32dvr(vin,resid,work)
     END IF
  ELSE IF(precon == 'block') THEN
     u0=1
     eig0=u0+maxblk*maxblk
     CALL rblvec(resid,vec,array(eig0),eig,array(u0),n,addvec)
  END IF
  IF(prnt) THEN
     title='new trial vectors iteration = '//itoc(iter)
     CALL prntrm(title,vec,n,addvec,n,addvec,iout)
  END IF
END SUBROUTINE genvec





