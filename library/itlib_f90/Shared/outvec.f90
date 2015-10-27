!deck outvec.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-06  Time: 08:29:56
 
!***begin prologue     outvec
!***date written       980420   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           write out vectors
!***author             schneider, barry (nsf)
!***source
!***
!***references

!***routines called
!***end prologue       outvec

SUBROUTINE outvec(vec,TYPE,n,nvc,prnt)

REAL*8, INTENT(IN OUT)                   :: vec(n,nvc)
CHARACTER (LEN=*), INTENT(IN OUT)        :: TYPE
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: nvc
INTEGER, INTENT(IN)                      :: prnt
IMPLICIT INTEGER (a-z)

CHARACTER (LEN=4) :: itoc
CHARACTER (LEN=80) :: title


COMMON/io/inp, iout

DO  i=1,nvc
  CALL iosys('write real "'//TYPE//itoc(i)//'" to rwf', n,vec(1,i),0,' ')
END DO
IF(prnt) THEN
  title='output vectors'
  CALL prntrm(title,vec,n,nvc,n,nvc,iout)
END IF
RETURN
END SUBROUTINE outvec

