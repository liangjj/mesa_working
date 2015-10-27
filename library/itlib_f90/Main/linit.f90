!deck linit.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-13  Time: 16:03:41
 
!***begin prologue     linit
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            initialize davidson
!***
!***description

!***references
!***routines called
!***end prologue       linit

SUBROUTINE linit(ht,hx,hy,hz,v,driver,vec,hvec,h,htmp,b,btmp,  &
    n,m,nt,nx,ny,nz,nbeg,nend,nout,spdim,maxvec)

REAL*8, INTENT(IN OUT)                   :: ht(nt,nt)
REAL*8, INTENT(IN OUT)                   :: hx(nx,nx)
REAL*8, INTENT(IN OUT)                   :: hy(ny,ny)
REAL*8, INTENT(IN OUT)                   :: hz(nz,nz)
REAL*8, INTENT(IN OUT)                   :: v(n,nt)
REAL*8, INTENT(IN OUT)                   :: driver(m)
REAL*8, INTENT(IN OUT)                   :: vec(m,maxvec)
REAL*8, INTENT(IN OUT)                   :: hvec(m,maxvec)
REAL*8, INTENT(IN OUT)                   :: h(maxvec,maxvec)
REAL*8, INTENT(IN OUT)                   :: htmp(maxvec,maxvec)
REAL*8, INTENT(IN OUT)                   :: b(maxvec)
REAL*8, INTENT(IN OUT)                   :: btmp(maxvec)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: m
INTEGER, INTENT(IN OUT)                  :: nt
INTEGER, INTENT(IN OUT)                  :: nx
INTEGER, INTENT(IN OUT)                  :: ny
INTEGER, INTENT(IN OUT)                  :: nz
INTEGER, INTENT(IN OUT)                  :: nbeg
INTEGER, INTENT(IN OUT)                  :: nend
INTEGER, INTENT(IN OUT)                  :: nout
INTEGER, INTENT(IN OUT)                  :: spdim
INTEGER, INTENT(IN OUT)                  :: maxvec
IMPLICIT INTEGER (a-z)






COMMON/io/inp, iout


!    initialize the effect of the hamiltonian on these vectors.

CALL hvdvr(ht,hx,hy,hz,v,vec(1,nbeg),hvec(1,nbeg),n,nt,nx,ny,nz, nout,spdim)

!     initialize the small hamiltonian matrix and right hand side.

CALL hinit(h,htmp,b,btmp,vec,hvec,driver,m,1,nend,maxvec)
RETURN
END SUBROUTINE linit



