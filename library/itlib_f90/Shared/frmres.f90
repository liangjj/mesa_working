!deck frmres.f
!***begin prologue     frmres
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           residual calculation
!***author             schneider, barry (nsf)
!***source
!***
!***references

!***routines called
!***end prologue       frmres
  SUBROUTINE frmres(rep,size,num2do,con,uncon,it,rtdone,code)
  USE io
  USE dvd_prnt
  USE dvd_global
  IMPLICIT NONE
  REAL*8                                 :: rep
  INTEGER                                :: size
  INTEGER                                :: num2do
  INTEGER                                :: con
  INTEGER                                :: uncon
  INTEGER                                :: it
  INTEGER                                :: rtdone
  INTEGER                                :: num2chk
  CHARACTER (LEN=*), DIMENSION(2)        :: code
  REAL*8                                 :: sdot, ERR, temp
  CHARACTER (LEN=16)                     :: STATUS
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=4)                      :: itoc
!        first transform the vectors to the new basis.

  CALL ebcxx(resid,vec,svec,n_dvd,size,size,n_dvd,n_dvd,maxvec)
  vec(1:n_dvd,1:size) = resid((1:n_dvd,1:size) 
!
!        second transform the effect of the hamiltonian on the basis
!        to the new basis.
  CALL ebcxx(resid,hvec,svec,n_dvd,size,size,n_dvd,n_dvd,maxvec)
  hvec(1:n_dvd,1:size) =  resid(1:n_dvd,1:size)
  IF(log_dvd(6)) THEN
     title='information for iteration = '//itoc(it)
     WRITE(iout,1) title
  END IF
!     form residuals for the desired targeted nroots.
  n2chk=MIN(size,num2do)
  DO  i=1,size
      resid(:,i) = hvec(:,i) - vec(:,i) * eigwrk(i)      
  END DO
  IF(log_dvd(7)) THEN
      title='transformed vectors iteration = '//itoc(it)
      CALL prntfm(title,vec,n_dvd,n2chk,n_dvd,maxvec,iout)
  END IF
  IF(log_dvd(8)) THEN
     title='hamiltonian on transformed vectors iteration = ' //itoc(it)
     CALL prntfm(title,hvec,n_dvd,n2chk,n_dvd,maxvec,iout)
  END IF
  IF(log_dvd(9)) THEN
     title='residuals iteration = '//itoc(it)
     CALL prntfm(title,resid,n_dvd,n2chk,n_dvd,maxvec,iout)
  END IF
!     re-constitute the small matrix
  b = 0.d0
  bwrk = 0.d0
  DO  i=1,size
      b(i,i) = eigwrk(i)
      bwrk(i,i) = eigwrk(i)
  END DO
!     check converged and unconverged roots
  uncon=0
  con=0
  maxerr=0.d0
  DO  i=1,n2chk
      ERR = SQRT (sdot(n_dvd,resid(1,i),1, resid(1,i),1) )
      temp=eigwrk(i) + rep
      maxerr=MAX(ERR,maxerr)
      IF(ERR <= cnverg) THEN
         STATUS='converged'
         con=con+1
         eig(i+rtdone)=eigwrk(i)
!           write out converged eigenpairs to rwf
         CALL iosys('write real "'//code(1)//itoc(rtdone+i)  &
                                 //'" to rwf',1,temp,0,' ')
         CALL iosys('write real "'//code(2)//itoc(rtdone+i)  &
                                 //'" to rwf',n,vec(1,i),0,' ')
         WRITE(iout,2) i, temp, ERR, STATUS
      ELSE
!           move eigenvalues and residuals of unconverged eigenpairs
         STATUS='unconverged'
         uncon=uncon+1
         resid(:,uncon) = resid(:,i)
         eigwrk(uncon) = eigwrk(i)
         WRITE(iout,2) i, temp, ERR, STATUS
      END IF
  END DO
1    FORMAT(/,5X,a80)
2    FORMAT(/,5X,'root            = ',i4,/,5X,  &
    'davidson energy = ',f15.8,/,5X, 'rms error       = ',f15.8,/,5X,  &
    'status          = ',a16)
END SUBROUTINE frmres

