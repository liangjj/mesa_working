!deck lindvd.f
!***begin prologue     lindvd
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           iterative eigenvalue
!***author             schneider, barry (nsf)
!***source
!***purpose            iterative linear system solver using
!***                   Davidson algorithm.
!***
!***description        solve the linear set of equations
!***                   ( E - H ) |Psi> = |B>
!***                   where |B> is a matrix.
!***references

!***routines called
!***end prologue       lindvd

SUBROUTINE lindvd(ibuf,hbuf,diag,energy,pvec,hpvec,  &
    rhs,h,htmp,b,btmp,resid,t,list,cnverg,thresh,  &
    n,nrhs,maxit,maxvec,lenbuf,incore,header, prnt)

INTEGER, INTENT(IN OUT)                  :: ibuf(2,lenbuf)
REAL*8, INTENT(IN OUT)                   :: hbuf(lenbuf)
REAL*8, INTENT(IN OUT)                   :: diag(n)
REAL*8, INTENT(IN OUT)                   :: energy
REAL*8, INTENT(IN OUT)                   :: pvec(n,maxvec)
REAL*8, INTENT(IN OUT)                   :: hpvec(n,maxvec)
REAL*8, INTENT(IN OUT)                   :: rhs(n,nrhs)
REAL*8, INTENT(IN OUT)                   :: h(maxvec,maxvec)
REAL*8, INTENT(IN OUT)                   :: htmp(maxvec,maxvec)
REAL*8, INTENT(IN OUT)                   :: b(maxvec,nrhs)
REAL*8, INTENT(IN OUT)                   :: btmp(maxvec,nrhs)
REAL*8, INTENT(IN OUT)                   :: resid(n,maxvec)
REAL*8, INTENT(IN OUT)                   :: t(n,*)
INTEGER, INTENT(OUT)                     :: list(nrhs)
REAL*8, INTENT(IN OUT)                   :: cnverg
REAL*8, INTENT(IN OUT)                   :: thresh
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: nrhs
INTEGER, INTENT(IN OUT)                  :: maxit
INTEGER, INTENT(IN)                      :: maxvec
INTEGER, INTENT(IN OUT)                  :: lenbuf
LOGICAL, INTENT(IN OUT)                  :: incore
CHARACTER (LEN=*), INTENT(IN OUT)        :: header(3)
LOGICAL, INTENT(IN)                      :: prnt(11)
IMPLICIT INTEGER (a-z)
  REAL*8                                 :: error
  REAL*8                                 :: sdot
  REAL*8                                 :: maxerr
  REAL*8                                 :: zero=0.d0
  REAL*8                                 :: one=1.d0

CHARACTER (LEN=5) :: itoc

CHARACTER (LEN=80) :: title

DO  i=1,nrhs
  list(i)=i
END DO
CALL iosys('read integer '//header(3)//' from hamiltonian', 1,ntot,0,' ')
WRITE(iout,1) nrhs, maxit, maxvec, cnverg
!-----------------------------------------------------------------------c
!                                                                       c
!                    Initialization Section                             c
!                                                                       c
!-----------------------------------------------------------------------c

!     initialize the set of input vectors as the
!     orthonormalized set of right hand side vectors

CALL copy(rhs,pvec,n*nrhs)
nbeg=1
nend=nrhs
CALL gschmt(pvec,thresh,n,1,nend,nout,.true.,.false.)
IF(nend == 0) THEN
  WRITE(iout,2)
  CALL lnkerr('quit linearsys')
ELSE
  nend=nout
  IF(prnt(2)) THEN
    title='initial trial vectors'
    CALL prntfm(title,pvec(1,nbeg),n,nout,n,nout,iout)
  END IF
END IF

!    initialize the effect of the hamiltonian on these vectors.

title='h on initial vectors'
Call packed_symmetric_matrix_on_vector((packed_matrix_d, ibuf, matrix_diagonal_d, vec_d(:,nbeg),    &
                                        h_vectors_d(:,nbeg), number)

CALL honv(ibuf,hbuf,diag,pvec(1,nbeg),hpvec(1,nbeg),  &
    n,nout,header,lenbuf,ntot,incore,title,prnt(3))

!     initialize the small hamiltonian matrix and right hand side.

CALL bsmall(h,htmp,b,btmp,pvec,hpvec,rhs,energy,nrhs,n,ndum,nend,  &
    maxvec,'initialize')
IF(prnt(4)) THEN
  title='initial small matrix'
  CALL prntfm(title,htmp,nend,nend,maxvec,maxvec,iout)
  title='initial small right hand side'
  CALL prntfm(title,btmp,nend,nrhs,maxvec,maxvec,iout)
END IF
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
newrhs=nrhs
totrhs=0
DO WHILE ( error > cnverg.AND.iter < maxit )
  iter = iter + 1
  WRITE(iout,4) iter, nend
!     Step 1:
  
!           get solutions of the small matrix.
  
!                btmp holds the initial matrix which is destroyed.
  
  
  CALL lsolve(htmp,btmp,resid,nend,newrhs,maxvec)
  IF(prnt(5)) THEN
    title='solutions of small matrix iteration = '//itoc(iter)
    CALL prntfm(title,btmp,nend,newrhs,maxvec,newrhs,iout)
  END IF
  
!     Step 2:
  
!           form the residuals and check for convergence.
!           t contains the all the solution vectors.
  
  
!           when lares is finished, the converged solutions are stored
!           on the disk and the unconverged residuals and solutions appear
!           first in the list.  the b and btmp matrices are updated to
!           reflect this and the current value of newrhs changed to the
!           new number.
  
  CALL lares(pvec,hpvec,b,btmp,rhs,energy,cnverg,resid,  &
      maxerr,t,list,n,nend,newrhs,con,uncon,maxvec, iter,prnt(5))
  newrhs=uncon
  totrhs=totrhs+con
  IF(totrhs == nrhs) THEN
    
!           all solutions are converged.  copy them in to rhs and quit
    
    RETURN
  ELSE
    
!           all solutions are not converged.  set the error to the largest
!           current error and continue the iteration sequence
    
    error=MIN(error,maxerr)
    
!           scale the residuals to get the next set of vectors to
!           be added to the krylov sequence.  leave them in resid
!           for the moment.
    
    CALL lavec(energy,diag,resid,n,uncon,iter,prnt(9))
    
!           orthogonalize the scaled residuals to the old vectors
!           to get an additional nout set of linearly independent
!           vectors.
    
    CALL abschm(pvec,resid,thresh,n,nend,uncon,nout, .true.,.false.)
    
!           orthonormalize these vectors and add them to the
!           existing set.
    
    CALL gschmt(resid,thresh,n,1,nout,nfinal,.true.,.false.)
    
!           how many new vectors could be added in principle
    
    numnew = maxvec - nend
    
!           how many will we add
    
    addvec = MIN(numnew,nfinal)
    WRITE(iout,5) nfinal, addvec
    
!           add them
    
    nbeg = nend + 1
    nend=nend+addvec
    IF(nend >= maxvec.AND.addvec == 0) THEN
      WRITE(iout,6)
      CALL lnkerr('vector space exceeded. quit')
    END IF
    CALL copy(resid,pvec(1,nbeg),n*addvec)
    IF(prnt(10)) THEN
      CALL tstovl(pvec,n,nend)
    END IF
    
!           we now have either a reconstituted set of vectors or
!           an augmented set of vectors.
    
!              operate with hamiltonian on these vectors
    
    title='h on initial vectors'
    CALL honv(ibuf,hbuf,diag,pvec(1,nbeg),hpvec(1,nbeg),  &
        n,addvec,header,lenbuf,ntot,incore, title,prnt(3))
    
!           update the curent small hamiltonian matrix and right hand sides.
    
    CALL bsmall(h,htmp,b,btmp,pvec,hpvec,rhs,energy,newrhs,  &
        n,nbeg,nend,maxvec,'fill')
  END IF
END DO
CALL iosys('write integer "size of davidson vector '//  &
    'space" to hamiltonian',1,nend,0,' ')
CALL iosys('write real "davidson vectors" to hamiltonian', nend*n,pvec,0,' ')
RETURN
1    FORMAT(/,1X,'davidson linear system solver',/,10X,  &
    'number of right hand sides    = ',i4,/,10X,  &
    'maximum number of iterations  = ',i4,/,10X,  &
    'maximum number of vectors     = ',i4,/,10X,  &
    'convergence criterion for rms = ',e15.8)
2    FORMAT(/,5X,'cannot even begin davidson calculation:',/,5X,  &
    'orthonormalization of initial vectors yields null' '  set')
3    FORMAT(/,5X,'beginning iterations:',/,5X, 'error initialized at ',e15.8)
4    FORMAT(/,5X,'iteration            = ',i5,  &
    /,5X,'size of vector space = ',i5)
5    FORMAT(/,5X,'number of linearly independent vectors after '  &
    'schmidt step = ',i4,/,5X,'number actually added ' '             = ',i4)
6    FORMAT(/,5X,'size of allowed vector space exceeded. quit')
END SUBROUTINE lindvd





