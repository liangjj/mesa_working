!deck dvd_err
!***begin prologue     dvd_err
!***date written       991013   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           error print and exit for davidson code.
!***author             schneider, barry (nsf)
!***source
!***purpose
!***
!***description
!***references

!***routines called
!***end prologue       dvd_err
  SUBROUTINE dvd_err
  USE io
  USE dvd_global
  IMPLICIT NONE
  IF(errno == 1) THEN
     WRITE(iout,1)
  ELSE IF(errno == 2) THEN
     WRITE(iout,2)
  ELSE IF(errno == 3) THEN
     WRITE(iout,3)
  END IF
  CALL lnkerr('quit davidson')
1    FORMAT(/,5X,'cannot even begin davidson calculation:',/,5X,  &
    'orthonormalization of initial vectors yields null' '  set')
2    FORMAT(/1X,'no more orthonormal vectors can be added')
3    FORMAT(/,5X,'cannot continue since maxvec is exceeded')
END SUBROUTINE dvd_err

