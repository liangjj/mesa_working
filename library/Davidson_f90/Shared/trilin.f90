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

SUBROUTINE trilin(vec,TYPE,n,ntrial,maxvec,prnt)

REAL*8, INTENT(IN OUT)                   :: vec(n,*)
CHARACTER (LEN=*), INTENT(IN OUT)        :: TYPE
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: ntrial
INTEGER, INTENT(IN)                      :: maxvec
INTEGER, INTENT(IN)                      :: prnt
IMPLICIT INTEGER (a-z)

CHARACTER (LEN=4) :: itoc
CHARACTER (LEN=80) :: title


COMMON/io/inp, iout

CALL iosys('create real trials on rwf',n*ntrial,0,0,' ')
start=0
DO WHILE(start < ntrial)
  END=MIN(start+maxvec,ntrial)
  start=start+1
  cnt=0
DO  i=start,END
cnt=cnt+1
CALL iosys('read real "'//TYPE//itoc(i)//'" from rwf', n,vec(1,cnt),0,' ')
END DO
CALL iosys('write real trials to rwf without rewinding', n*cnt,vec,0,' ')
CALL iosys('rewind trials on rwf read-and-write',0,0,0,' ')
IF(prnt) THEN
  title='input vectors'
  CALL prntrm(title,vec,n,cnt,n,cnt,iout)
END IF
start=END
END DO
RETURN
END SUBROUTINE trilin

