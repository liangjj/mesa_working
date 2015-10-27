!deck prepit
!***begin prologue     prepit
!***date written       010828   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           iterative eigenvalue
!***author             schneider, barry (nsf)
!***source
!***purpose            prepare for eigdvr routine
!***references

!***routines called
!***end prologue       prepit

SUBROUTINE prepit(pham,px,pv,dim,n,nd,nroot,ntrial,nattim,  &
    cnverg,thresh,precon,maxit,maxvec,nblck,dvdprt)

INTEGER*8, INTENT(IN)                    :: pham(*)
INTEGER*8, INTENT(IN OUT)                :: px(*)
INTEGER*8, INTENT(IN)                    :: pv(dim+1)
INTEGER, INTENT(IN)                      :: dim
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: nd(dim+1)
INTEGER, INTENT(IN)                      :: nroot
INTEGER, INTENT(IN)                      :: ntrial
INTEGER, INTENT(IN OUT)                  :: nattim
REAL*8, INTENT(IN OUT)                   :: cnverg
REAL*8, INTENT(IN OUT)                   :: thresh
CHARACTER (LEN=*), INTENT(IN)            :: precon
INTEGER, INTENT(IN OUT)                  :: maxit
INTEGER, INTENT(IN)                      :: maxvec
INTEGER, INTENT(IN)                      :: nblck
LOGICAL, INTENT(IN OUT)                  :: dvdprt(12)
IMPLICIT INTEGER (a-z)
CHARACTER (LEN=8) :: key
REAL*8 hx, hy, hz, vpot, trial, h, dvd

#ifdef decpointer

INTEGER*8 phx, phy, phz, pre
#END IF decpointer
#ifdef sgipointer

INTEGER*4 phx, phy, phz, pre
#END IF sgipointer
LOGICAL :: zeroit, prnt



COMMON/io/inp, iout
pointer (ptrial,trial(1)), (ptrial,itrial(1))
pointer (pdvd,dvd(1))
pointer (phx,hx(1))
pointer (phy,hy(1))
pointer (phz,hz(1))
pointer (pvpot,vpot(1))
pointer (pre,h(1))

!     calculate the perturbation potential

key='$vpert'
zeroit=.true.
CALL vpert(key,vword,zeroit,prnt)
pvpot=pv(dim+1)

phx=pham(1)
phy=pham(2)
phz=pham(3)
vec=1
eigv=vec+n*ntrial
need=wpadti(eigv+ntrial)

IF(precon /= 'block') THEN
  
!        prepare for preconditioning and trial vectors based on the separable
!        hamiltonian ( hx + hy + hz )
  
  CALL seppre(pre,hx,hy,hz,nd(1),nd(2),nd(3),n,dim)
  IF(dim >= 2) THEN
    ind=need
    need=need+n*dim
  END IF
  CALL getmem(need,ptrial,ngot,'trial',0)
  IF(dim == 2) THEN
    CALL setnd2(itrial(ind),nd(1),nd(2),n)
  END IF
  IF(dim == 3) THEN
    CALL setnd3(itrial(ind),nd(1),nd(2),nd(3),n)
  END IF
  u1=1
  eig1=u1+nd(1)*nd(1)
  u2=eig1+nd(1)
  eig2=u2+nd(2)*nd(2)
  u3=eig2+nd(2)
  eig3=u3+nd(3)*nd(3)
  CALL vtrsep(h(u1),h(u2),h(u3),h(eig1),h(eig2),h(eig3),  &
      trial(vec),trial(eigv),itrial(ind), nd(1),nd(2),nd(3),dim,n,ntrial,.false.)
  
ELSE IF(precon == 'block') THEN
  
!        prepare for preconditioning based on block diagonalization
!        of full hamiltonian.
  
  CALL blkpre(pre,hx,hy,hz,vpot,nd(1),nd(2),nd(3),n,nblck,dim)
  u=1
  eig=u+nblck*nblck
  ind=need
  need=need+n
  CALL getmem(need,ptrial,ngot,'trial',0)
  CALL vtrblk(h(eig),h(u),trial(eigv),trial(vec),itrial(ind),  &
      n,ntrial,nblck,.false.)
END IF

!     the trial vectors are stored, we can get rid of the memory.

CALL getmem(-ngot,ptrial,idum,'trial',idum)

!     now set up the memory for the davidson subroutine.

diag=1
eig=diag+n
vec=eig+nroot
hvec=vec+n*maxvec
resid=hvec+n*maxvec
b=resid+n*maxvec
bwrk=b+maxvec*maxvec
eigwrk=bwrk+maxvec*maxvec
work=eigwrk+maxvec
svec=work+MAX(5*maxvec,n*maxvec)
need=wpadti(svec+maxvec*maxvec)
CALL getmem(need,pdvd,nword,'dvd',0)

!     store the full diagonal element and then zero
!     the diagonal elements of the one-dimensional matrices.

CALL diagnl(dvd(diag),hx,hy,hz,vpot,n,nd(1),nd(2),nd(3),dim)

!     the potential is no longer needed, its on the diagonal.

CALL getmem(-vword,pvpot,ndum,'vint',ndum)

CALL dvrdvd(hx,hy,hz,dvd(diag),dvd(eig),dvd(vec),dvd(hvec),  &
    dvd(resid),dvd(b),dvd(bwrk),dvd(eigwrk),dvd(work),  &
    dvd(svec),cnverg,thresh,dim,nd(1),nd(2),nd(3),  &
    n,nroot,ntrial,nattim,maxit,maxvec,nblck, dvdprt,precon,pre)
CALL getmem(-nword,pdvd,idum,'dvd',idum)
RETURN
END SUBROUTINE prepit


