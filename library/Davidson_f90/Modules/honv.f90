!deck h_on_v.f
!***begin prologue     h_on_v
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           
!***author             schneider, barry (nsf)
!***source
!***purpose            multiply the hamiltonian matrix on a vector.
!***references

!***routines called
!***end prologue       h_on_v

SUBROUTINE honv(ibuf,hbuf,diag,veco,vecn,n,nvc,headr,  &
    lenbuf,nel,incore,title,prnt)

INTEGER, INTENT(IN)                      :: ibuf(2,lenbuf)
REAL*8, INTENT(IN)                       :: hbuf(lenbuf)
REAL*8, INTENT(IN)                       :: diag(n)
REAL*8, INTENT(IN)                       :: veco(n,nvc)
REAL*8, INTENT(OUT)                      :: vecn(n,nvc)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: nvc
CHARACTER (LEN=*), INTENT(IN OUT)        :: headr(3)
INTEGER, INTENT(IN)                      :: lenbuf
INTEGER, INTENT(IN)                      :: nel
LOGICAL, INTENT(IN)                      :: incore
CHARACTER (LEN=*), INTENT(IN OUT)        :: title
LOGICAL, INTENT(IN)                      :: prnt
IMPLICIT INTEGER (a-z)
CALL rzero(vecn,n*nvc)
IF(nel /= 0) THEN
  IF(incore) THEN
    DO  i=1,nel
      ii=ibuf(1,i)
      jj=ibuf(2,i)
      DO  j=1,nvc
        vecn(ii,j) = vecn(ii,j) + hbuf(i)*veco(jj,j)
        vecn(jj,j) = vecn(jj,j) + hbuf(i)*veco(ii,j)
      END DO
    END DO
  ELSE
    CALL iosys('rewind all on hamiltonian',0,0,0,' ')
    trips=nel/lenbuf
    left=nel-trips*lenbuf
    DO  nt=1,trips
      CALL iosys('read integer '//headr(2)//' from '//  &
          'hamiltonian without rewinding', 2*lenbuf,ibuf,0,' ')
      CALL iosys('read integer '//headr(2)//' from '//  &
          'hamiltonian without rewinding', wptoin(lenbuf),hbuf,0,' ')
      DO  i=1,lenbuf
        ii=ibuf(1,i)
        jj=ibuf(2,i)
        DO  j=1,nvc
          vecn(ii,j) = vecn(ii,j) + hbuf(i)*veco(jj,j)
          vecn(jj,j) = vecn(jj,j) + hbuf(i)*veco(ii,j)
        END DO
      END DO
    END DO
    IF(left /= 0) THEN
      CALL iosys('read integer '//headr(2)//' from '//  &
          'hamiltonian without rewinding', 2*left,ibuf,0,' ')
      CALL iosys('read integer '//headr(2)//' from '//  &
          'hamiltonian without rewinding', wptoin(left),hbuf,0,' ')
      DO  i=1,left
        ii=ibuf(1,i)
        jj=ibuf(2,i)
        DO  j=1,nvc
          vecn(ii,j) = vecn(ii,j) + hbuf(i)*veco(jj,j)
          vecn(jj,j) = vecn(jj,j) + hbuf(i)*veco(ii,j)
        END DO
      END DO
    END IF
  END IF
END IF
DO  i=1,n
  DO  j=1,nvc
    vecn(i,j) = vecn(i,j) + diag(i)*veco(i,j)
  END DO
END DO
IF(prnt) THEN
  CALL prntrm(title,vecn,n,nvc,n,nvc,iout)
END IF
RETURN
END SUBROUTINE honv
