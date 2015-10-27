!*deck dvd_dat
!***begin prologue     dvd_dat
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           
!***author             schneider, barry (nsf)
!***source             
!***purpose            data entry for dvr davidson routine.
!***                   
!***references         
!
!***routines called    
!***end prologue       dvd_dat
SUBROUTINE dvd_dat
  USE io
  USE dvd_prnt
  USE dvd_global
  REAL*8                           :: fpkey
  CHARACTER (LEN=80)               :: chrkey
  LOGICAL                          :: dollar, logkey
!
  if( dollar('$dvr_dvd',card,cpass,inp) ) then 
      cnverg=fpkey(card,'convergence',1.d-08,' ')
      thresh=fpkey(card,'overlap-tolerance',1.d-08,' ')
      precon=chrkey(card,'preconditioner','none',' ')
      if(precon.ne.'none') then
         nblck=intkey(card,'maximum-size-of-preconditioning-block', &
                      200,' ')                      
      endif
      nroot=intkey(card,'number-of-roots',n,' ')
      nroot=min(nroot,n)
      maxvec=intkey(card,'maximum-number-of-vectors',2*nroot,' ')
      maxvec=min(maxvec,n)
      maxit=intkey(card,'maximum-number-of-iterations',maxvec,' ')
      nattim=intkey(card,'number-of-roots-at-a-time',nroot,' ')
      ntrial=intkey(card,'number-of-trial-vectors',nroot,' ')
      ntrial = max(ntrial,nroot)
      ntrial=min(ntrial,maxvec)
      prn_dvd_loc = 'print=davidson='//prn_dvd
      prn_dvd_loc(12) = chrkey(card,'print=davidson=',prn_dvd(12),' ')
      IF (prn_dvd_loc(12) == 'all') THEN
          log_dvd = .true.
      ELSE
          CALL setlog(log_dvd,prn_dvd_loc,card,11)
      END IF
  endif
  write(iout,1) nroot, nattim, ntrial, &
                thresh, cnverg, maxit, maxvec, precon, nblck
1 format(/,15x,'iterative diagonalization information',/,/,5x,      &
               'number of roots                    = ',i3,/,5x,     &
               'number of roots at a time          = ',i3,/,5x,     &
               'number of trials                   = ',i3,/,5x,     &  
               'overlap tolerance                  = ',e15.8,/,5x,  &
               'convergence criterion              = ',e15.8,/,5x,  &
               'maximum number of iterations       = ',i6,/,5x,     &
               'maximum number of vectors          = ',i6,/,5x,     &       
               'preconditioning                    = ',a24,/,5x,    &
               'block size                         = ',i5)           
END SUBROUTINE dvd_dat




