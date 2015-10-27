!deck trilin.f
!***begin prologue     trilin
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           read in and write out trial vectors
!***author             schneider, barry (nsf)
!***source
!***
!***references

!***routines called
!***end prologue       trilin

  SUBROUTINE trilin(TYPE)
  USE io
  USE dvd_prnt
  USE dvd_global
  IMPLICIT NONE
  CHARACTER (LEN=*), INTENT(IN)        :: TYPE
  INTEGER                              :: start, end, cnt, i
  CHARACTER (LEN=4)                    :: itoc
!
  CALL iosys('create real trials on rwf',n_dvd*ntrial,0,0,' ')
  start=0
  DO WHILE(start < ntrial)
     END=MIN(start+maxvec,ntrial)
     start=start+1
     cnt=0
     DO  i=start,END
         cnt=cnt+1
         CALL iosys('read real "'//TYPE//itoc(i)//'" from rwf', &
                     n_dvd,vec(1,cnt),0,' ')
     END DO
     CALL iosys('write real trials to rwf without rewinding',   &
                 n_dvd*cnt,vec,0,' ')
     CALL iosys('rewind trials on rwf read-and-write',0,0,0,' ')
     IF(log_dvd(1)) THEN
        title='input vectors'
        CALL prntrm(title,vec,n_dvd,cnt,n_dvd,cnt,iout)
     END IF
     start=END
  END DO
END SUBROUTINE trilin

