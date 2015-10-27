!deck lindvr.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-13  Time: 16:03:38
 
!***begin prologue     lindvr
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           iterative eigenvalue
!***author             schneider, barry (nsf)
!***source
!***purpose            iterative linear system solver specialized
!***                   to time-DVR using the Davidson algorithm.
!***
!***description        solve the linear set of equations
!***                        A |Psi> = |B>
!***                   where A and |B> are matrices.
!***
!***                   m = nz*ny*nx*nt*nc*2
!***                   n = nz*ny*nx*nt
!***references

!***routines called
!***end prologue       lindvr

SUBROUTINE lindvr(hx,hy,hz,ht,eigx,eigy,eigz,ux,uy,uz,  &
    v,pvec,hpvec,h,htmp,b,btmp,rhs,trials,resid,  &
    soln,t1,t2,lufac,luind,scale,cnverg,thresh,eps,  &
    nx,ny,nz,nt,nc,dim,m,n,nrhs,ntrial,maxit, maxvec,prnt,precon)

REAL*8, INTENT(IN OUT)                   :: hx(nx,nx)
REAL*8, INTENT(IN OUT)                   :: hy(ny,ny)
REAL*8, INTENT(IN OUT)                   :: hz(nz,nz)
REAL*8, INTENT(IN OUT)                   :: ht(nt,nt)
REAL*8, INTENT(IN OUT)                   :: eigx(nx)
REAL*8, INTENT(IN OUT)                   :: eigy(ny)
REAL*8, INTENT(IN OUT)                   :: eigz(nz)
REAL*8, INTENT(IN OUT)                   :: ux(nx,nx)
REAL*8, INTENT(IN OUT)                   :: uy(ny,ny)
REAL*8, INTENT(IN OUT)                   :: uz(nz,nz)
REAL*8, INTENT(IN OUT)                   :: v(*)
REAL*8, INTENT(IN OUT)                   :: pvec(m,maxvec)
REAL*8, INTENT(IN OUT)                   :: hpvec(m,maxvec)
REAL*8, INTENT(IN OUT)                   :: h(maxvec,maxvec)
REAL*8, INTENT(IN OUT)                   :: htmp(maxvec,maxvec)
REAL*8, INTENT(IN OUT)                   :: b(maxvec,nrhs)
REAL*8, INTENT(IN OUT)                   :: btmp(maxvec,nrhs)
REAL*8, INTENT(IN OUT)                   :: rhs(m,nrhs)
REAL*8, INTENT(IN OUT)                   :: trials(m,ntrial)
REAL*8, INTENT(IN OUT)                   :: resid(m,maxvec)
REAL*8, INTENT(IN OUT)                   :: soln(m,nrhs)
REAL*8, INTENT(IN OUT)                   :: t1(*)
REAL*8, INTENT(IN OUT)                   :: t2(*)
COMPLEX*16, INTENT(IN OUT)               :: lufac(*)
INTEGER, INTENT(IN OUT)                  :: luind(*)
REAL*8, INTENT(OUT)                      :: scale
REAL*8, INTENT(IN OUT)                   :: cnverg
REAL*8, INTENT(IN OUT)                   :: thresh
REAL*8, INTENT(IN)                       :: eps
INTEGER, INTENT(IN OUT)                  :: nx
INTEGER, INTENT(IN OUT)                  :: ny
INTEGER, INTENT(IN OUT)                  :: nz
INTEGER, INTENT(IN OUT)                  :: nt
INTEGER, INTENT(IN OUT)                  :: nc
INTEGER, INTENT(IN)                      :: dim
INTEGER, INTENT(IN OUT)                  :: m
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: nrhs
INTEGER, INTENT(IN)                      :: ntrial
INTEGER, INTENT(IN OUT)                  :: maxit
INTEGER, INTENT(IN)                      :: maxvec
LOGICAL, INTENT(IN)                      :: prnt(11)
CHARACTER (LEN=*), INTENT(IN)            :: precon
IMPLICIT INTEGER (a-z)


REAL*8  error, sdot
REAL*8 maxerr
REAL*4 time(20), delta(20), secnds


CHARACTER (LEN=5) :: itoc
CHARACTER (LEN=1) :: it
CHARACTER (LEN=80) :: title











COMMON/io/inp, iout
WRITE(iout,1) m, nrhs, maxit, maxvec, ntrial, cnverg, thresh, eps
WRITE(iout,2) scale
!-----------------------------------------------------------------------c
!                                                                       c
!                    Initialization Section                             c
!                                                                       c
!-----------------------------------------------------------------------c

!     initialize the set of input vectors as the

IF(prnt(2)) THEN
  title='input trial vectors'
  CALL prntrm(title,trials,m,ntrial,m,ntrial,iout)
  title='right hand side'
  CALL prntrm(title,rhs,m,nrhs,m,nrhs,iout)
  title='initial solution'
  CALL prntrm(title,soln,m,nrhs,m,nrhs,iout)
END IF
CALL copy(trials,pvec,m*ntrial)
time(1)=secnds(0.0)
nbeg=1
nend=ntrial
CALL gschmt(pvec,thresh,m,1,nend,nout,.true.,.false.)
time(2)=secnds(0.0)
delta(1)=time(2)-time(1)
IF(nend == 0) THEN
  CALL dvderr(1)
ELSE
  nend=nout
  IF(prnt(2)) THEN
    title='initial trial vectors'
    CALL prntfm(title,pvec(1,nbeg),m,nout,m,nout,iout)
  END IF
END IF
!    initialize the effect of the hamiltonian on these vectors.

title='h on initial vectors'
CALL htonv(hx,hy,hz,ht,v,pvec(1,nbeg),hpvec(1,nbeg),n,nc,  &
    nx,ny,nz,nt,nout,dim)
time(3)=secnds(0.0)
delta(2)=time(3)-time(2)

!     initialize the small hamiltonian matrix and right hand side.

CALL hinit(h,htmp,b,btmp,pvec,hpvec,rhs,m,nrhs,nend,maxvec)
time(4)=secnds(0.0)
delta(3)=time(4)-time(3)
IF(prnt(4)) THEN
  title='initial small matrix'
  CALL prntfm(title,htmp,nend,nend,maxvec,maxvec,iout)
  title='initial small right hand side'
  CALL prntfm(title,btmp,nend,nrhs,maxvec,maxvec,iout)
END IF
WRITE(iout,7) (delta(i),i=1,3)
!----------------------------------------------------------------------c
!                                                                      c
!                    Iteration Sequence                                c
!                                                                      c
!     iteration is continued until all of the solution are converged   c
!     or if convergence is not achieved some preset maximum number of  c
!     iterations are performed.                                        c
!                                                                      c
!----------------------------------------------------------------------c
iter=0
error=1.d+10
WRITE(iout,3) error
DO WHILE ( error > cnverg.AND.iter < maxit )
  time(5)=secnds(0.0)
  iter = iter + 1
  WRITE(iout,4) iter, nend
!     Step 1:
  
!           get solutions of the small matrix.
  
!                btmp holds the initial matrix which is destroyed.
!                note that resid is used as temporary storage in vscale.
  
  title='iteration = '//itoc(iter)//' solving linear system '  &
      //'of dimension N = '//itoc(nend)
  WRITE(iout,5) title
  
  CALL lslv(htmp,btmp,resid,nend,nrhs,maxvec)
  CALL newsol(soln,pvec,btmp,m,nend,nrhs,maxvec)
  IF(prnt(5)) THEN
    title='solutions of small matrix iteration = '//itoc(iter)
    CALL prntfm(title,btmp,nend,nrhs,maxvec,nrhs,iout)
  END IF
  time(6)=secnds(0.0)
  delta(4)=time(6)-time(5)
  
!     Step 2:
  
!           form the residuals and check for convergence.
!           t1 contains the transformed vectors and t2 the transformed
!           hamiltonian on vectors.
  
  
!        unconverged residuals are in resid and the corresponding eigenvalues
!        are in etmp.
  
  CALL lres(pvec,hpvec,btmp,rhs,scale,cnverg,resid,  &
      maxerr,soln,t1,m,nend,nrhs,con,uncon, maxvec,iter,prnt(5))
  time(7)=secnds(0.0)
  delta(5)=time(7)-time(6)
  IF(maxerr <= eps) THEN
    WRITE(iout,*) 'writing guess vectors to ham'
    CALL iosys('write real "guess vectors" to ham', nrhs*m,t1,0,' ')
  END IF
  IF(con == nrhs) THEN
!            call copy(t1,soln,m*nrhs)
    CALL copy(soln,rhs,m*nrhs)
    
!           all solutions are converged.  copy them in to rhs and quit
    
    title='final solution'
    CALL prntrm(title,rhs,m,nrhs,m,nrhs,iout)
    RETURN
  ELSE
    
!           all solutions are not converged.  set the error to the largest
!           current error and continue the iteration sequence or
!           quit if maxvec is exceeded.
    
    error=MIN(error,maxerr)
!           how many new vectors could be added in principle
    
    numnew = maxvec - nend
    
!           how many will we add
    
    addvec = MIN(numnew,uncon)
    WRITE(iout,6) addvec
    chkno = nend + addvec
    IF(chkno <= maxvec) THEN
      
!              maximum number of vectors is still within the allowed
!              limits.  add vectors to the set from the unconverged
!              residuals and put them after the current vectors.
      
      nbeg = nend + 1
      IF(precon == 'none') THEN
        IF(dim == 1) THEN
          CALL v1t2h0(resid,t1,ux,nx,nt,nc,addvec)
          CALL v1te(t1,eigx,nx,nt,nc,addvec)
          CALL v1t2dvr(t1,pvec(1,nbeg),ux,nx,nt,nc,addvec)
        ELSE IF(dim == 2) THEN
          CALL v2t2h0(resid,t1,ux,uy,nx,ny,nt,nc,addvec)
          CALL v2te(t1,eigx,eigy,nx,ny,nt,nc,addvec)
          CALL v2t2dvr(t1,pvec(1,nbeg),ux,uy,nx,ny,nt, nc,addvec)
        ELSE IF(dim == 3) THEN
          CALL v3t2h0(resid,t1,ux,uy,uz,nx,ny,nz,nt, nc,addvec)
          CALL v3te(t1,eigx,eigy,eigz,nx,ny,nz,nt,nc,addvec)
          CALL v3t2dvr(t1,pvec(1,nbeg),ux,uy,uz,nx,ny,nz, nt,nc,addvec)
        END IF
      ELSE
        CALL cblslv(resid,pvec(1,nbeg),t1,lufac,luind,n,nc, addvec)
      END IF
      time(8)=secnds(0.0)
      delta(6)=time(8)-time(7)
      nend=nend+addvec
      
!              orthonormalize the new trials to the old vectors
!              to get an additional nout vectors.
      
      CALL gschmt(pvec,thresh,m,nbeg,nend,nout,.true.,.false.)
      time(9)=secnds(0.0)
      delta(7)=time(9)-time(8)
      
      IF(nout == 0) THEN
        
!                 no more vectors write out unconverged results.
        
!                  call dvderr(2)
        RETURN
      END IF
      nend=nbeg+nout-1
!               if(prnt(10)) then
!                  call tstovl(pvec,m,nend)
!               endif
      
!           we now have either a reconstituted set of vectors or
!           an augmented set of vectors.
      
!              operate with hamiltonian on these vectors
      
      title='h on initial vectors'
      CALL htonv(hx,hy,hz,ht,v,pvec(1,nbeg),hpvec(1,nbeg),  &
          n,nc,nx,ny,nz,nt,nout,dim)
      time(10)=secnds(0.0)
      delta(8)=time(10)-time(9)
      
!           update the curent small hamiltonian matrix and right hand sides.
      
      CALL hupdat(h,htmp,b,btmp,pvec,hpvec,rhs,m,nrhs,nbeg, nend,maxvec)
      time(11)=secnds(0.0)
      delta(9)=time(11)-time(10)
    ELSE
      CALL dvderr(3)
    END IF
  END IF
  WRITE(iout,8) iter, (delta(i),i=4,9)
END DO
CALL iosys('write integer "size of davidson vector '//  &
    'space" to ham',1,nend,0,' ')
CALL iosys('write real "davidson vectors" to ham', nend*m,pvec,0,' ')
RETURN
1    FORMAT(/,1X,'davidson linear system solver using preconditioning',  &
    /,10X, 'number of equations                 = ',i4,/,10X,  &
    'number of right hand sides          = ',i4,/,10X,  &
    'maximum number of iterations        = ',i4,/,10X,  &
    'maximum number of vectors           = ',i4,/,10X,  &
    'number of initial trial vectors     = ',i4,/,10X,  &
    'convergence criterion for rms       = ',e15.8,/,10X,  &
    'overlap threshold criterion         = ',e15.8,/,10X,  &
    'restart criterion                   = ',e15.8)
2    FORMAT(/,5X,'hamiltonian scale factor = ',e15.8)
3    FORMAT(/,5X,'beginning davidson iterations:',/,5X,  &
    'error initialized at ',e15.8)
4    FORMAT(/,1X,'beginning next davidson cycle',/,5X,  &
    'iteration = ',i4,1X,'size of vector space = ',i5)
5    FORMAT(/,5X,a80)
6    FORMAT(/,5X,'number of vectors actually added = ',i4)
7    FORMAT(/,1X,'time for schmidt orthogonalization of initial '  &
    'vectors = ',e15.8, /,1X,'time for H on initial vectors          '  &
    '               = ',e15.8,  &
    /,1X,'time for initialization of small  hamiltonian       ' '  = ',e15.8)
8    FORMAT(/,1X,'timing at iteration = ',i4,  &
    /,1X,'time for solution of small set of linear ' 'equations = ',e15.8,  &
    /,1X,'time for residual calculation                      ' '= ',e15.8,  &
    /,1X,'time for preconditioning                         ' '  = ',e15.8,  &
    /,1X,'time for schmidt orthogonalization              ' '   = ',e15.8,  &
    /,1X,'time for H on vectors                           ' '   = ',e15.8,  &
    /,1X,'time to update small matrix                     ' '   = ',e15.8)
END SUBROUTINE lindvr

