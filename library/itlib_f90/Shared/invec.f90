!deck invec.f
!***begin prologue     invec
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           read in vectors
!***author             schneider, barry (nsf)
!***source
!***
!***references

!***routines called
!***end prologue       invec

  SUBROUTINE invec(TYPE,nvc)
  USE io
  USE dvd_prnt
  USE dvd_global
  IMPLICIT NONE
  CHARACTER (LEN=*)                      :: TYPE
  INTEGER                                :: nvc
  CHARACTER (LEN=4)                      :: itoc
  CHARACTER (LEN=80)                     :: title
  DO  i=1,nvc
      CALL iosys('read real "'//TYPE//itoc(i)//'" from rwf',   &
                  n_dvd,resid(1,i),0,' ')
  END DO
  IF(log_dvd(1)) THEN
     title='input vectors'
     CALL prntrm(title,vec,n_dvd,nvc,n_dvd,nvc,iout)
  END IF
END SUBROUTINE invec

