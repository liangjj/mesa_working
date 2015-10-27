!deck lres.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-13  Time: 16:03:50
 
!***begin prologue     lres
!***date written       980420   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           residual calculation
!***author             schneider, barry (nsf)
!***source
!***
!***references

!***routines called
!***end prologue       lres

SUBROUTINE lres(vec,hvec,coef,rhs,scale,cnverg,resid,  &
    maxerr,soln,t1,n,m,nrhs,ncon,addvec,maxvec, it,prnt)

REAL*8, INTENT(IN OUT)                   :: vec(n,*)
REAL*8, INTENT(IN OUT)                   :: hvec(n,*)
REAL*8, INTENT(IN OUT)                   :: coef(maxvec,nrhs)
REAL*8, INTENT(IN)                       :: rhs(n,nrhs)
REAL*8, INTENT(IN)                       :: scale
REAL*8, INTENT(IN)                       :: cnverg
REAL*8, INTENT(OUT)                      :: resid(n,*)
REAL*8, INTENT(OUT)                      :: maxerr
REAL*8, INTENT(IN OUT)                   :: soln(n,*)
REAL*8, INTENT(IN)                       :: t1(n,*)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(IN)                      :: nrhs
INTEGER, INTENT(OUT)                     :: ncon
INTEGER, INTENT(OUT)                     :: addvec
INTEGER, INTENT(IN)                      :: maxvec
INTEGER, INTENT(IN OUT)                  :: it
LOGICAL, INTENT(IN)                      :: prnt(4)
IMPLICIT INTEGER (a-z)

REAL*8  sdot, ERR
CHARACTER (LEN=16) :: STATUS
CHARACTER (LEN=80) :: title
CHARACTER (LEN=4) :: itoc



COMMON/io/inp, iout

!        calculate the effect of the hamiltonian on the solution vectors
!        and put them in t1

CALL ebcxx(t1,hvec,coef,n,m,nrhs,n,n,maxvec)
title='information for iteration = '//itoc(it)
WRITE(iout,1) title
DO  i=1,nrhs
  DO  j=1,n
    resid(j,i) = rhs(j,i) - t1(j,i)
  END DO
END DO
IF(prnt(2)) THEN
  title='residuals iteration = '//itoc(it)
  CALL prntfm(title,resid,n,nrhs,n,maxvec,iout)
END IF
IF(prnt(3)) THEN
  title='solution vectors iteration = '//itoc(it)
  CALL prntfm(title,soln,n,nrhs,n,maxvec,iout)
END IF
IF(prnt(4)) THEN
  title='hamiltonian on solution vectors iteration = ' //itoc(it)
  CALL prntfm(title,t1,n,nrhs,n,maxvec,iout)
END IF
addvec=0
ncon=0
maxerr=0.d0
DO  i=1,nrhs
  ERR = scale*SQRT (sdot(n,resid(1,i),1, resid(1,i),1) )
  maxerr=MAX(ERR,maxerr)
  IF(ERR <= cnverg) THEN
    STATUS='converged'
    ncon=ncon+1
    CALL iosys('write real "solution for right hand side = '  &
        //itoc(i)//'" to ham',n,soln(1,i),0,' ')
    WRITE(iout,2) i, ERR, STATUS
  ELSE
    STATUS='unconverged'
    addvec=addvec+1
    CALL copy(resid(1,i),resid(1,addvec),n)
    WRITE(iout,2) i, ERR, STATUS
  END IF
END DO
test=m+addvec
IF(test > maxvec) THEN
  addvec=maxvec-m
  WRITE(iout,3) addvec
END IF
RETURN
1    FORMAT(/,5X,a80)
2    FORMAT(/,5X,'solution        = ',i4,/,5X,  &
    'rms error       = ',e15.8,/,5X, 'status          = ',a16)
3    FORMAT(/,5X,'number of added vectors will exceed maxvec',/,5X,  &
    'number of vectors actually added = ',i4)
END SUBROUTINE lres













