!deck timpts.f
!***begin prologue     timpts
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            time-dependent potential
!***
!***description        automat time grid
!***
!***references

!***routines called
!***end prologue       timpts

  SUBROUTINE timpts(t,npts,ntreg)
  USE dvr_global,    ONLY   : input, output
  IMPLICIT NONE
  INTEGER                               :: ntreg
  REAL*8, DIMENSION(ntreg+1)            :: t(ntreg+1)
  INTEGER, DIMENSION(ntreg)             :: npts(ntreg)
  REAL*8  delt
  LOGICAL :: posinp
  CHARACTER (LEN=1600) :: card
  INTEGER :: i, j
  IF ( posinp('$pts',card) ) THEN
       READ(input,*) t(1), delt, npts(1)
  ELSE
       CALL lnkerr('not found')
  END IF
  DO  i=2,ntreg+1
    t(i)=t(i-1)+delt
  END DO
  DO  i=2,ntreg
    npts(i)=npts(i-1)
  END DO
  RETURN
END SUBROUTINE timpts


