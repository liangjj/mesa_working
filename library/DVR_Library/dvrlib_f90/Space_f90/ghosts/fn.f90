!deck fn.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:43:29
 
!***begin prologue     fn
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose
!***
!***references

!***routines called
!***end prologue       fn

SUBROUTINE fn(pt,f,n,ftyp,typ)

REAL*8, INTENT(IN)                       :: pt(n)
REAL*8, INTENT(OUT)                      :: f(n)
INTEGER, INTENT(IN)                      :: n
CHARACTER (LEN=*), INTENT(IN)            :: ftyp
CHARACTER (LEN=*), INTENT(IN)            :: typ
IMPLICIT INTEGER (a-z)



COMMON/io/inp, iout

IF(ftyp == 'negative-exponential') THEN
  IF(typ == 'function') THEN
    DO  i=1,n
      f(i)=EXP(-pt(i))
    END DO
  ELSE IF(typ == 'first-derivative') THEN
    DO  i=1,n
      f(i)=- EXP(-pt(i))
    END DO
  ELSE IF(typ == 'second-derivative') THEN
    DO  i=1,n
      f(i)=EXP(-pt(i))
    END DO
  END IF
ELSE IF(ftyp == 'positive-exponential') THEN
  IF(typ == 'function') THEN
    DO  i=1,n
      f(i)=EXP(pt(i))
    END DO
  ELSE IF(typ == 'first-derivative') THEN
    DO  i=1,n
      f(i)= EXP(pt(i))
    END DO
  ELSE IF(typ == 'second-derivative') THEN
    DO  i=1,n
      f(i)=EXP(pt(i))
    END DO
  END IF
ELSE IF(ftyp == 'sine') THEN
  IF(typ == 'function') THEN
    DO  i=1,n
      f(i)=SIN(pt(i))
    END DO
  ELSE IF(typ == 'first-derivative') THEN
    DO  i=1,n
      f(i)= COS(pt(i))
    END DO
  ELSE IF(typ == 'second-derivative') THEN
    DO  i=1,n
      f(i)= - SIN(pt(i))
    END DO
  END IF
ELSE IF(ftyp == 'cosine') THEN
  IF(typ == 'function') THEN
    DO  i=1,n
      f(i)=COS(pt(i))
    END DO
  ELSE IF(typ == 'first-derivative') THEN
    DO  i=1,n
      f(i)= - SIN(pt(i))
    END DO
  ELSE IF(typ == 'second-derivative') THEN
    DO  i=1,n
      f(i)= - COS(pt(i))
    END DO
  END IF
ELSE IF(ftyp == 'x**2') THEN
  IF(typ == 'function') THEN
    DO  i=1,n
      f(i)=pt(i)*pt(i)
    END DO
  ELSE IF(typ == 'first-derivative') THEN
    DO  i=1,n
      f(i)= 2.d0*pt(i)
    END DO
  ELSE IF(typ == 'second-derivative') THEN
    DO  i=1,n
      f(i)= 2.d0
    END DO
  END IF
ELSE IF(ftyp == 'x') THEN
  IF(typ == 'function') THEN
    DO  i=1,n
      f(i)=pt(i)
    END DO
  ELSE IF(typ == 'first-derivative') THEN
    DO  i=1,n
      f(i)= 1.d0
    END DO
  ELSE IF(typ == 'second-derivative') THEN
    DO  i=1,n
      f(i)= 0.d0
    END DO
  END IF
END IF
RETURN
END SUBROUTINE fn



