!deck @(#)guess2.f 1.1  11/30/90
 
! Code converted using TO_F90 by Alan Miller
! Date: 2008-03-08  Time: 10:54:36

SUBROUTINE guess2(ptgues,nwks,ibuf,rbuf,lnbuf,  &
    hguess,nwksg,diag,t1,eigvec,eigval, ops,ntotal,repcor,prtflg,dagtyp,title)

!***begin prologue     guess2
!***date written       870813   (yymmdd)
!***revision date      yymmdd   (yymmdd)

!***keywords           ci, guess
!***author             saxe, paul (lanl)
!***source             @(#)guess2.f 1.1   11/30/90

!***purpose            to form and diagonalize a portion of the h
!     matrix to obtain guesses at ci vectors.

!***description

!***references

!***routines called    (none)

!***end prologue       guess2


INTEGER, INTENT(OUT)                     :: ptgues(nwks)
INTEGER, INTENT(IN)                      :: nwks
INTEGER, INTENT(IN OUT)                  :: ibuf(2,lnbuf)
REAL*8, INTENT(IN)                       :: rbuf(lnbuf)
INTEGER, INTENT(IN)                      :: lnbuf
REAL*8, INTENT(OUT)                      :: hguess(nwksg,nwksg)
INTEGER, INTENT(IN)                      :: nwksg
REAL*8, INTENT(IN)                       :: diag(nwks)
REAL*8, INTENT(IN OUT)                   :: t1(nwks)
REAL*8, INTENT(IN)                       :: eigvec(nwksg,nwksg)
REAL*8, INTENT(IN OUT)                   :: eigval(nwksg)
CHARACTER (LEN=*), INTENT(IN OUT)        :: ops
INTEGER, INTENT(IN)                      :: ntotal
REAL*8, INTENT(IN OUT)                   :: repcor
CHARACTER (LEN=*), INTENT(IN)            :: prtflg
LOGICAL, INTENT(IN OUT)                  :: dagtyp
CHARACTER (LEN=*), INTENT(IN OUT)        :: title(3)
IMPLICIT INTEGER (a-z)



CHARACTER (LEN=80) :: htit




INTEGER :: nnpg






REAL*8 eguess
REAL*8 t


CHARACTER (LEN=4) :: itoc
INTEGER :: intkey
LOGICAL :: logkey


COMMON /io/ inp,iout

!     ----- fill an array pointing to guess numbers, or zero if unused

CALL iosys('read real '//title(3)//' from hamiltonian', nwks,t1,0,' ')
CALL copy(t1,diag,nwks)
CALL izero(ptgues,nwks)

DO  i=1,nwksg
  eguess=1.0D+30
  DO  j=1,nwks
    IF(t1(j) >= eguess) GO TO 10
    refwlk=j
    eguess=t1(j)
    10         CONTINUE
  END DO
  t1(refwlk)=1.0D+30
  ptgues(refwlk)=i
END DO

!     ----- read through the large hamiltonian and extract the small

CALL rzero(hguess,nwksg*nwksg)

npass=ntotal/lnbuf
left = ntotal - npass*lnbuf
CALL iosys('rewind all on hamiltonian',0,0,0,' ')

DO  pass=1,npass
  CALL iosys('read integer '//title(1)//' from hamiltonian'//  &
      ' without rewinding',2*lnbuf,ibuf,0,' ')
  CALL iosys('read integer '//title(1)//' from hamiltonian'//  &
      ' without rewinding',wptoin(lnbuf),rbuf,0,' ')
  
  DO  n=1,lnbuf
    i=ptgues(ibuf(1,n))
    IF (i <= 0) CYCLE
    j=ptgues(ibuf(2,n))
    IF (j <= 0) CYCLE
    ii=MAX(i,j)
    jj=MIN(i,j)
    hguess(ii,jj)=hguess(ii,jj)+rbuf(n)
  END DO
END DO

IF(left /= 0) THEN
  CALL iosys('read integer '//title(1)//' from hamiltonian'//  &
      ' without rewinding',2*left,ibuf,0,' ')
  CALL iosys('read integer '//title(1)//' from hamiltonian'//  &
      ' without rewinding',wptoin(left),rbuf,0,' ')
  DO  n=1,left
    i=ptgues(ibuf(1,n))
    IF (i <= 0) CYCLE
    j=ptgues(ibuf(2,n))
    IF (j <= 0) CYCLE
    ii=MAX(i,j)
    jj=MIN(i,j)
    hguess(ii,jj)=hguess(ii,jj)+rbuf(n)
  END DO
END IF

!     ----- add in the diagonals -----

DO  j=1,nwks
  i=ptgues(j)
  IF (i <= 0) CYCLE
  hguess(i,i)=diag(j)

















END DO

DO  i=1,nwksg
  DO  j=1,i
    hguess(j,i)=hguess(i,j)
  END DO
END DO
IF (logkey(ops,'print=ci=guess-matrix',.false.,' ')) THEN
  htit='guess hamiltonian'
  CALL prntrm(htit,hguess,nwksg,nwksg,nwksg,nwksg,iout)
END IF

!        ----- diagonalize -----

CALL dsyev('v','l',nwksg,hguess,nwksg,eigval,t1,5*nwksg,info)

IF (info /= 0) THEN
  CALL lnkerr('error in diagonalization')
END IF

nroots=intkey(ops,'ci=nroots',1,' ')
nroots=MIN(nwksg,nroots)

IF (prtflg /= 'minimum') THEN
  WRITE (iout,1) nwksg
END IF

DO  root=1,nwksg
  
!        ----- expand up the vector -----
  
  CALL rzero(t1,nwks)
  
  t=0.0D+00
  DO  i=1,nwks
    IF (ptgues(i) > 0) THEN
      t1(i)=eigvec(ptgues(i),root)
      IF (ABS(t1(i)) > t) THEN
        t=ABS(t1(i))
        ref=i
      END IF
    END IF
  END DO
  
  IF (prtflg /= 'minimum') THEN
    WRITE (iout,2) root,ref,diag(ref)+repcor, eigval(root)+repcor, t1(ref)
  END IF
  
  CALL iosys('write real "guess:'//itoc(root)//'" to guess', nwks,t1,0,' ')
END DO
CALL iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')


RETURN
1    FORMAT (/,t15,'results from guess matrix of size',i5,/,  &
    ' guess   reference   diagonal energy  ' 'guess energy        c(0)')
2    FORMAT (1X,i3,i10,4X,g20.9,g16.9,f8.4)
END SUBROUTINE guess2
