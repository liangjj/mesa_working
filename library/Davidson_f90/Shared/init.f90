! \documentclass{article}
 
! Code converted using TO_F90 by Alan Miller
! Date: 2008-03-16  Time: 12:02:59
 
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{INIT: Driver for Initialization of Davidson Vectors/Matrix}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle

!deck init.f
!***begin prologue     init
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           initialize, davidson
!***purpose            driver for initialization of davidson vectors
!***                   and matrices.
!***description
!***author             schneider, barry (nsf)
!***source
!***
!***references

!***routines called
!***end prologue       init

SUBROUTINE init(hx,hy,hz,diag,vec,hvec,b,bwrk,scr,thresh,  &
    drctv,dim,n,nin,begin,size,nx,ny,nz,rtdone, code,maxvec,orth,prnt)

REAL*8, INTENT(IN OUT)                   :: hx(nx,nx)
REAL*8, INTENT(IN OUT)                   :: hy(ny,ny)
REAL*8, INTENT(IN OUT)                   :: hz(nz,nz)
REAL*8, INTENT(IN OUT)                   :: diag(n)
REAL*8, INTENT(IN OUT)                   :: vec(n,*)
REAL*8, INTENT(IN OUT)                   :: hvec(n,*)
REAL*8, INTENT(IN OUT)                   :: b(maxvec,*)
REAL*8, INTENT(IN OUT)                   :: bwrk(maxvec,*)
REAL*8, INTENT(IN OUT)                   :: scr(*)
REAL*8, INTENT(IN OUT)                   :: thresh
CHARACTER (LEN=*), INTENT(IN OUT)        :: drctv
INTEGER, INTENT(IN OUT)                  :: dim
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: nin
INTEGER, INTENT(IN)                      :: begin
INTEGER, INTENT(OUT)                     :: size
INTEGER, INTENT(IN OUT)                  :: nx
INTEGER, INTENT(IN OUT)                  :: ny
INTEGER, INTENT(IN OUT)                  :: nz
INTEGER, INTENT(IN)                      :: rtdone
CHARACTER (LEN=*), INTENT(IN OUT)        :: code
INTEGER, INTENT(IN OUT)                  :: maxvec
LOGICAL, INTENT(IN)                      :: orth
LOGICAL, INTENT(IN OUT)                  :: prnt(3)
IMPLICIT INTEGER (a-z)


CHARACTER (LEN=80) :: title
CHARACTER (LEN=4) :: itoc






COMMON/io/inp, iout

!     get the additional set of vectors by orthogonalizing the trials
!     to the set of roots which have been converged.  this will ensure
!     that the subspace is orthogonal to previously converged vectors.

nout = nin
IF(orth) THEN
  IF(rtdone /= 0) THEN
    CALL invec(scr,code,n,rtdone,prnt(1))
    CALL abschm(scr,vec(1,begin),thresh,n,rtdone,nin, nout,.true.,.false.)
  END IF
END IF
END = begin + nout - 1
IF(nout /= 0) THEN
  CALL gschmt(vec,thresh,n,begin,END,nout,.true.,prnt(2))
END IF
IF(nout == 0) THEN
  CALL lnkerr('quit davidson. no more trial vectors '// 'possible')
END IF
size = begin + nout - 1

!     initialize the effect of the hamiltonian on these vectors.

title='h on initial vectors'
CALL hvdvr(hx,hy,hz,diag,vec(1,begin),hvec(1,begin),  &
    n,nx,ny,nz,nout,dim,prnt(3))

!     initialize the small hamiltonian matrix.

CALL hsmall(b,bwrk,vec,hvec,n,begin,size,maxvec,drctv,.false.)
RETURN
END SUBROUTINE init

