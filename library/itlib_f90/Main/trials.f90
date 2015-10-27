!deck trials.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-13  Time: 16:04:32
 
!***begin prologue     trials
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           iterative eigenvalue
!***author             schneider, barry (nsf)
!***source
!***purpose            iterative linear system solver specialized
!***                   to time-DVR using the Davidson algorithm.
!***
!***description        calculate an orthonormal set of trial vectors
!***
!***                   m = nz*ny*nx*nt*2
!***                   nxyzt = nz*ny*nx*nt
!***                   nxyz = nz*ny*nx
!***references

!***routines called
!***end prologue       trials

SUBROUTINE trials(trial,rhs,thresh,nx,ny,nz,nt,ntrial, m,n,dim,prn)

REAL*8, INTENT(IN OUT)                   :: trial(m,n)
REAL*8, INTENT(IN OUT)                   :: rhs(m,*)
REAL*8, INTENT(IN OUT)                   :: thresh
INTEGER, INTENT(IN OUT)                  :: nx
INTEGER, INTENT(IN OUT)                  :: ny
INTEGER, INTENT(IN OUT)                  :: nz
INTEGER, INTENT(IN OUT)                  :: nt
INTEGER, INTENT(OUT)                     :: ntrial
INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: dim
LOGICAL, INTENT(IN)                      :: prn
IMPLICIT INTEGER (a-z)


CHARACTER (LEN=800) :: card
CHARACTER (LEN=3) :: ans
CHARACTER (LEN=80) :: cpass, chrkey, typetr
LOGICAL :: dollar

COMMON /io/ inp, iout

IF ( dollar('$trials',card,cpass,inp) ) THEN
  typetr=chrkey(card,'type-of-trial-vectors','unit',' ')
  ntrial=intkey(card,'number-of-trial-vectors',m,' ')
  ntrial=MIN(ntrial,n)
END IF
IF(typetr == 'unit') THEN
  CALL rzero(trial,m*ntrial)
  add=MAX(m/2,1)
  DO  i=1,ntrial
    trial(i,i)=1.d0
    trial(i+add,i)=1.d0
  END DO
ELSE IF(typetr == 'guess-vectors') THEN
  CALL iosys('does "guess vectors" exist on ham',0,0,0,ans)
  IF(ans == 'yes') THEN
    CALL iosys('read real "guess vectors" from ham',m*ntrial, trial,0,' ')
  ELSE
    CALL lnkerr('trial vectors do not exist')
  END IF
ELSE
  CALL lnkerr('no option for other type trial vectors yet')
END IF
CALL gschmt(trial,thresh,m,1,ntrial,nout,.true.,.false.)
ntrial=nout
IF(prn) THEN
  cpass='trial vectors for '//typetr//' trials'
  CALL prntrm(cpass,trial,m,ntrial,m,ntrial,iout)
END IF
WRITE(iout,1) ntrial, typetr
1    FORMAT(/,10X,'number of trial vectors = ',i4,/,10X,  &
    'type of trial vectors   = ',a24)
RETURN
END SUBROUTINE trials

