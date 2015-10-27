! \documentclass{article}
 
! Code converted using TO_F90 by Alan Miller
! Date: 2008-03-16  Time: 12:02:39
 
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{FRMRES: Form Residuals and Transform to Optimal Vectors}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck frmres.f
!***begin prologue     frmres
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           residual calculation
!***author             schneider, barry (nsf)
!***source
!***
!***references

!***routines called
!***end prologue       frmres

SUBROUTINE frmres(vec,hvec,resid,trmat,eig,eigwrk,b,bwrk,rep,  &
    cnverg,maxerr,n,m,nroot,con,uncon, maxvec,it,prnt,point,code)

REAL*8, INTENT(IN)                       :: vec(n,*)
REAL*8, INTENT(IN)                       :: hvec(n,*)
REAL*8, INTENT(OUT)                      :: resid(n,*)
REAL*8, INTENT(IN OUT)                   :: trmat(maxvec,*)
REAL*8, INTENT(OUT)                      :: eig(*)
REAL*8, INTENT(IN OUT)                   :: eigwrk(*)
REAL*8, INTENT(OUT)                      :: b(maxvec,*)
REAL*8, INTENT(OUT)                      :: bwrk(maxvec,*)
REAL*8, INTENT(IN)                       :: rep
REAL*8, INTENT(IN)                       :: cnverg
REAL*8, INTENT(OUT)                      :: maxerr
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(IN OUT)                  :: nroot
INTEGER, INTENT(OUT)                     :: con
INTEGER, INTENT(OUT)                     :: uncon
INTEGER, INTENT(IN OUT)                  :: maxvec
INTEGER, INTENT(IN OUT)                  :: it
LOGICAL, INTENT(IN)                      :: prnt(4)
INTEGER, INTENT(IN OUT)                  :: point
CHARACTER (LEN=*), INTENT(IN OUT)        :: code(2)
IMPLICIT INTEGER (a-z)

REAL*8  sdot, ERR, temp
CHARACTER (LEN=16) :: STATUS
CHARACTER (LEN=80) :: title
CHARACTER (LEN=4) :: itoc





COMMON/io/inp, iout

n2chk=MIN(m,nroot)

!        first transform the vectors to the new basis.

CALL ebcxx(resid,vec,trmat,n,m,m,n,n,maxvec)
CALL copy(resid,vec,n*m)

!        second transform the effect of the hamiltonian on the basis
!        to the new basis.

CALL ebcxx(resid,hvec,trmat,n,m,m,n,n,maxvec)
CALL copy(resid,hvec,n*m)
IF(prnt(1)) THEN
  title='information for iteration = '//itoc(it)
  WRITE(iout,1) title
END IF

!     form residuals for the desired targeted nroots.

DO  i=1,m
  DO  j=1,n
    resid(j,i) = hvec(j,i) - eigwrk(i)*vec(j,i)
  END DO
END DO
IF(prnt(2)) THEN
  title='transformed vectors iteration = '//itoc(it)
  CALL prntfm(title,vec,n,n2chk,n,maxvec,iout)
END IF
IF(prnt(3)) THEN
  title='hamiltonian on transformed vectors iteration = ' //itoc(it)
  CALL prntfm(title,hvec,n,n2chk,n,maxvec,iout)
END IF
IF(prnt(4)) THEN
  title='residuals iteration = '//itoc(it)
  CALL prntfm(title,resid,n,n2chk,n,maxvec,iout)
END IF

!     re-constitute the small matrix

CALL rzero(b,maxvec*maxvec)
CALL rzero(bwrk,maxvec*maxvec)
DO  i=1,m
  b(i,i) = eigwrk(i)
  bwrk(i,i) = eigwrk(i)
END DO

!     check converged and unconverged roots

uncon=0
con=0
maxerr=0.d0
DO  i=1,n2chk
  ERR = SQRT (sdot(n,resid(1,i),1, resid(1,i),1) )
  temp=eigwrk(i) + rep
  maxerr=MAX(ERR,maxerr)
  IF(ERR <= cnverg) THEN
    STATUS='converged'
    con=con+1
    eig(i+point)=eigwrk(i)
    
!           write out converged eigenpairs to rwf
    
    CALL iosys('write real "'//code(1)//itoc(point+i)  &
        //'" to rwf',1,temp,0,' ')
    CALL iosys('write real "'//code(2)//itoc(point+i)  &
        //'" to rwf',n,vec(1,i),0,' ')
    WRITE(iout,2) i, temp, ERR, STATUS
  ELSE
    
!           move eigenvalues and residuals of unconverged eigenpairs
    
    STATUS='unconverged'
    uncon=uncon+1
    CALL copy(resid(1,i),resid(1,uncon),n)
    eigwrk(uncon) = eigwrk(i)
    WRITE(iout,2) i, temp, ERR, STATUS
  END IF
END DO
RETURN
1    FORMAT(/,5X,a80)
2    FORMAT(/,5X,'root            = ',i4,/,5X,  &
    'davidson energy = ',f15.8,/,5X, 'rms error       = ',f15.8,/,5X,  &
    'status          = ',a16)
END SUBROUTINE frmres

