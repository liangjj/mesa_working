!deck dvddat.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2008-03-01  Time: 16:31:00
 
!***begin prologue     dvddat
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            data entry for davidson routine.
!***
!***references

!***routines called
!***end prologue       dvdddat

SUBROUTINE dvddat(ops,nwks,nwksg,nroots,nattim,cnverg,thresh,  &
    mxiter,mxdvd,nguess,left,prdvd,dvdall)

CHARACTER (LEN=*), INTENT(IN OUT)        :: ops
INTEGER, INTENT(IN OUT)                  :: nwks
INTEGER, INTENT(OUT)                     :: nwksg
INTEGER, INTENT(IN)                      :: nroots
INTEGER, INTENT(OUT)                     :: nattim
REAL*8, INTENT(OUT)                      :: cnverg
REAL*8, INTENT(OUT)                      :: thresh
INTEGER, INTENT(OUT)                     :: mxiter
INTEGER, INTENT(OUT)                     :: mxdvd
INTEGER, INTENT(OUT)                     :: nguess
INTEGER, INTENT(IN OUT)                  :: left
LOGICAL, INTENT(OUT)                     :: prdvd(11)
LOGICAL, INTENT(OUT)                     :: dvdall
IMPLICIT INTEGER (a-z)
REAL*8  fpkey

LOGICAL :: logkey

COMMON/io/inp, iout

thresh=fpkey(ops,'ci=tolerance',1.0D-10,' ')
cnverg=fpkey(ops,'ci=convergence',1.0D-08,' ')
nattim=intkey(ops,'ci=nroots-at-a-time',nroots,' ')
nwksg=intkey(ops,'ci=guess-size',0,' ')
nguess=intkey(ops,'ci=initital-vector-size',0,' ')


!         ----- core allocation for guess routine -----

IF (nwksg <= 0) THEN
  DO  nwksg=10,500,10
    IF (wpadti(nwksg*nwksg+nwksg+MAX(5*nwksg,nwks)) + nwks > left) GO TO 20
  END DO
  nwksg=500
  GO TO 30
  20           CONTINUE
  nwksg=nwksg-10
  30           CONTINUE
END IF

nwksg=MIN(nwks,nwksg)
nguess=MIN(nguess,nwksg)
nguess=MAX(nguess,nroots)
IF(logkey(ops,'mcscf',.false.,' ')) THEN
  mxiter=intkey(ops,'mcscf=ci=iterations',20,' ')
ELSE
  mxiter=intkey(ops,'ci=iterations',15,' ')
  mxdvd=intkey(ops,'ci=davidson-vectors',15,' ')
  mxdvd=MAX(mxdvd,3*nguess)
  mxdvd=MIN(mxdvd,nwks)
END IF
IF(mxiter < 2*nattim) THEN
  mxiter=2*nattim
END IF
mxiter=MIN(mxiter,nwks)
prdvd(1)=logkey(card,'print=ci=trials',.false.,' ')
prdvd(2)=logkey(card,'print=ci=vectors',.false.,' ')
prdvd(3)=logkey(card,'print=ci=h-on-vectors',.false.,' ')
prdvd(4)=logkey(card,'print=ci=hamiltonian',.false.,' ')
prdvd(5)=logkey(card,'print=ci=iteration-information',.false.,' ')
prdvd(6)=logkey(card,'print=ci=residuals',.false.,' ')
prdvd(7)=logkey(card,'print=ci=transformed-vectors',.false.,' ')
prdvd(8)=logkey(card,'print=ci=transformed-h-on-vectors', .false.,' ')
prdvd(9)=logkey(card,'print=ci=new-trial-vectors',.false.,' ')
prdvd(10)=logkey(card,'print=ci=overlaps',.false.,' ')
dvdall=logkey(card,'print=ci=all',.false.,' ')
prdvd(11)=.false.
IF(dvdall) THEN
  DO  i=1,11
    prdvd(i)=.true.
  END DO
END IF
WRITE(iout,1) nroots, nattim, nwks, nwksg, nguess, mxdvd, mxiter,  &
    thresh, cnverg
1    FORMAT(/,15X,'iterative diagonalization information',/,/,5X,  &
    'number of roots                     = ',i4,/,5X,  &
    'number of roots at a time           = ',i4,/,5X,  &
    'size of matrix                      = ',i8,/,5X,  &
    'size of guess matrix                = ',i7,/,5X,  &
    'size of initial vectors             = ',i7,/,5X,  &
    'maximum number of davidson vectors  = ',i7,/,5X,  &
    'max. number of iterations           = ',i7,/,5X,  &
    'overlap tolerance                   = ',e15.8,/,5X,  &
    'convergence criterion               = ',e15.8)

RETURN
END SUBROUTINE dvddat




