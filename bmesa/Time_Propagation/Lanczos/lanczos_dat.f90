!deck lanczos_dat
!**begin prologue     lanczos_dat
!**date written       960718   (yymmdd)
!**revision date               (yymmdd)
!**keywords           time, dvr, Lanczos, propagation
!**
!**author             schneider, b. i.(nsf)
!**source             Lanczos
!**purpose            set up the parameters specifying the maximum
!**                   number of Lanczos vectors, convergence criterion and
!***                  overlap tolerance for vectors kept.
!**references
!**routines called    iosys, util and mdutil
!**end prologue       lanczos_dat
  SUBROUTINE lanczos_dat
  USE lanczos_global
  USE lanczos_prnt
  IMPLICIT NONE
  LOGICAL                       :: dollar, logkey
  CHARACTER (LEN=80)            :: chrkey
  REAL*8                        :: fpkey
  INTEGER                       :: intkey 
  IF ( dollar('$lanczos_dat',card,cpass,inp) ) then
      pr_prp='print=lanczos='//pr_prp
      pr_prp(7)=chrkey(card,'print=lanczos=',pr_prp(7),' ')
      IF(pr_prp(7) == 'all') THEN
         CALL setprn(log_prp,6)
      ELSE
         CALL setlog(log_prp,pr_prp,card,6)
      END IF
      cnverg=fpkey(card,'convergence',1.d-08,' ')
      maxit=intkey(card,'maximum-number-of-iterations',n3d,' ')
      maxit=MIN(maxit,n3d)
      WRITE(iout,1) cnverg, maxit
  ELSE
      write(iout,2)
      stop
  END IF
1 FORMAT(/,15X,'lanczos diagonalization information',/,&
         /,5X,'convergence criterion              = ',e15.8, &
         /,5X,'maximum number of iterations       = ',i6)
2 FORMAT(/,5x,'error in diagonalization card section')
END SUBROUTINE lanczos_dat




