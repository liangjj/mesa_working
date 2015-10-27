!deck plymat.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:45:44
 
!***begin prologue     plymat
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            lobatto matrix elements.
!***
!***references

!***routines called
!***end prologue       plymat

SUBROUTINE plymat(q,wt,p,dp,ddp,ov,fder,sder,fdfd,bloch, n,coord,prn)

REAL*8, INTENT(IN OUT)                   :: q(n)
REAL*8, INTENT(IN OUT)                   :: wt(n)
REAL*8, INTENT(IN OUT)                   :: p(n,n)
REAL*8, INTENT(IN OUT)                   :: dp(n,n)
REAL*8, INTENT(IN OUT)                   :: ddp(n,n)
REAL*8, INTENT(IN OUT)                   :: ov(n,n)
REAL*8, INTENT(IN OUT)                   :: fder(n,n)
REAL*8, INTENT(IN OUT)                   :: sder(n,n)
REAL*8, INTENT(IN OUT)                   :: fdfd(n,n)
REAL*8, INTENT(IN OUT)                   :: bloch(n,n)
INTEGER, INTENT(IN OUT)                  :: n
CHARACTER (LEN=*), INTENT(IN)            :: coord
LOGICAL, INTENT(IN)                      :: prn
IMPLICIT INTEGER (a-z)



CHARACTER (LEN=80) :: title

CHARACTER (LEN=16) :: typarg

COMMON/io/inp, iout

typarg='linear'
IF(coord == 'cylindrical') THEN
  typarg='quadratic'
END IF
CALL fafb(p,p,wt,ov,n)
CALL fadfb(p,dp,q,wt,fder,typarg,n)
CALL faddfb(p,dp,ddp,q,wt,sder,typarg,n)
CALL dfadfb(dp,dp,q,wt,fdfd,typarg,n)
CALL tbloch(p,dp,q,bloch,typarg,n)
IF(prn) THEN
  title='overlap matrix'
  CALL prntrm(title,ov,n,n,n,n,iout)
  title='psi*dpsi matrix'
  CALL prntrm(title,fder,n,n,n,n,iout)
  title='psi*ddpsi matrix'
  CALL prntrm(title,sder,n,n,n,n,iout)
  title='dpsi*dpsi matrix'
  CALL prntrm(title,fdfd,n,n,n,n,iout)
  title='bloch matrix'
  CALL prntrm(title,bloch,n,n,n,n,iout)
  CALL maket(p,dp,ddp,q,wt,sder,bloch,fdfd,typarg,n)
  title='symmetrized kinetic energy matrix'
  CALL prntrm(title,fdfd,n,n,n,n,iout)
END IF
RETURN
END SUBROUTINE plymat






