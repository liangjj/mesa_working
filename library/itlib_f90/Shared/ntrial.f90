!deck ntrial.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-06  Time: 08:29:46
 
!***begin prologue     ntrial
!***date written       980420   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           vectors, rearrange
!***author             schneider, barry (nsf)
!***purpose            make a new set of trial vectors for
!***                   the next set of roots.
!***description
!                      vec     = vectors holding the trials
!                      n       = vector length
!                      size    = size of small matrix diagonalized in current
!                                pass
!                      con     = number of converged vectors on the
!                                current pass
!                      nin     = number of trial vectors used
!                      trial   = number of trial vectors available
!                      maxvec  = maximum number of allowable vectors

!***source             itlib
!***
!***references

!***routines called
!***end prologue       ntrial

SUBROUTINE ntrial(vec,n,size,con,nin,trial,maxvec)

REAL*8, INTENT(IN OUT)                   :: vec(n,*)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: size
INTEGER, INTENT(IN)                      :: con
INTEGER, INTENT(IN OUT)                  :: nin
INTEGER, INTENT(IN OUT)                  :: trial
INTEGER, INTENT(IN OUT)                  :: maxvec
IMPLICIT INTEGER (a-z)

CHARACTER (LEN=4) :: itoc

COMMON/io/inp, iout

!     there are,

navail = size - con

!     vectors from the diagonalization of the small matrix which are
!     rayleigh-ritz approximations to the next set of roots.

n2cpy=MIN(navail,trial,maxvec)
IF(n2cpy > 0) THEN
  CALL copy(vec(1,con+1),vec,n*n2cpy)
END IF

!     if there is still room left, add the unused trial vectors, provided
!     there are any left.

nout=n2cpy
toadd=trial- nin
IF(toadd > 0) THEN
  nout=MIN(n2cpy+toadd,trial,maxvec)
  add=nout-n2cpy
  nout=n2cpy+add
  IF(add > n2cpy) THEN
    point=nin+1
    DO  i=n2cpy+1,nout
      CALL iosys('read real "trial:'//itoc(point) //'" from rwf',  &
          n,vec(1,i),0,' ')
      point=point+1
    END DO
  END IF
END IF
WRITE(iout,1) n2cpy, add, nout
IF(nout /= 0) THEN
  DO  i=1,nout
    CALL iosys('write real "trial:'//itoc(i) //'" to rwf',  &
        n,vec(1,i),0,' ')
  END DO
ELSE
  WRITE(iout,2)
  CALL lnkerr('quit:no trials')
END IF

!     update values of trial and nin

trial=nout
nin=nout
RETURN
1    FORMAT(/,1X,'creating a new set of trial vectors',  &
    /,1X,'number of unconverged roots copied = ',i6,  &
    /,1X,'number of unused trials copied     = ',i6,  &
    /,1X,'total number of new trials         = ',i6)
2    FORMAT(/,1X,'no trial vectors are left:quit')
END SUBROUTINE ntrial


