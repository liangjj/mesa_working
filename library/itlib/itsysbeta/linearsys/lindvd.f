*deck lindvd.f
c***begin prologue     lindvd
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative eigenvalue
c***author             schneider, barry (nsf)
c***source             
c***purpose            iterative linear system solver using 
c***                   Davidson algorithm.
c***                   
c***description        solve the linear set of equations
c***                   ( E - H ) |Psi> = |B>
c***                   where |B> is a matrix.
c***references         
c
c***routines called    
c***end prologue       lindvd
      subroutine lindvd(ibuf,hbuf,diag,energy,pvec,hpvec,
     1                  rhs,h,htmp,b,btmp,resid,t,list,cnverg,thresh,
     2                  n,nrhs,maxit,maxvec,lenbuf,ntot,prnt)
      implicit integer (a-z)
      real*8 hbuf, diag, energy, pvec, hpvec, rhs
      real*8 h, htmp, b, btmp
      real*8 resid, t, error, cnverg, thresh, tmp
      real*8 maxerr, zero, one, nrzero
      logical prnt, incore
      character*5 itoc
      character*80 title
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n)
      dimension pvec(n,maxvec), hpvec(n,maxvec), rhs(n,nrhs)
      dimension h(maxvec,maxvec), htmp(maxvec,maxvec)
      dimension b(maxvec,nrhs), btmp(maxvec,nrhs)
      dimension resid(n,maxvec), t(n,*), list(nrhs)
      dimension prnt(11)
      common/io/inp, iout
      data zero, one/ 0.d0, 1.d0/
      data nrzero / 1.0d-06 /
      data incore/.false./
      data hamno/90/
      open(unit=70,file='gmres',err=99,form='unformatted',
     1     status='unknown')    
      do 10 i=1,nrhs
         list(i)=i
 10   continue   
      lwork=5*maxvec
c
c     if ntot is small enough, the non-zero elements and indices will be
c     assumed to reside in the buffers and not be read from disk over and
c     over again. 
c
      if(ntot.le.lenbuf) then
         incore=.true.
	 write(iout,*) 'Matrix Vector Multiply in Core'
      endif      
      write(iout,1) n, nrhs, maxit, maxvec, cnverg  
c-----------------------------------------------------------------------c
c                                                                       c
c                    Initialization Section                             c
c                                                                       c
c-----------------------------------------------------------------------c
c
c     initialize the set of input vectors as the 
c     orthonormalized set of right hand side vectors
c
      call copy(rhs,pvec,n*nrhs)
      nbeg=1
      nend=nrhs
      call gschmt(pvec,thresh,n,1,nend,nout,.true.,.false.)
      if(nend.eq.0) then
         write(iout,2)
         stop
      else
         nend=nout
         if(prnt(2)) then
            title='initial trial vectors'
            call prntfm(title,pvec(1,nbeg),n,nout,n,nout,iout)
         endif
      endif
c
c    initialize the effect of the hamiltonian on these vectors.
c     
      title='h on initial vectors'
      call honv(ibuf,hbuf,diag,pvec(1,nbeg),hpvec(1,nbeg),
     1          n,nout,lenbuf,ntot,incore,0,hamno,prnt(3))
c
c     initialize the small hamiltonian matrix and right hand side.
c
      call bsmall(h,htmp,b,btmp,pvec,hpvec,rhs,energy,nrhs,n,ndum,nend,
     1            maxvec,'initialize')       
      if(prnt(4)) then
         title='initial small matrix'
         call prntrm(title,htmp,nend,nend,maxvec,maxvec,iout)
         title='initial small right hand side'
         call prntrm(title,btmp,nend,nrhs,maxvec,maxvec,iout)
      endif
c----------------------------------------------------------------------c
c                                                                      c
c                    Iteration Sequence                                c
c                                                                      c
c     iteration is continued until all of the solution are converged   c 
c     or if convergence is not achieved some preset maximum number of  c
c     iterations are performed.                                        c
c                                                                      c
c----------------------------------------------------------------------c
      iter=0
      error=1.d+10
      write(iout,3) error
      newrhs=nrhs
      totrhs=0
      do while ( error.gt.cnverg.and.iter.lt.maxit )
         iter = iter + 1
         write(iout,4) iter, nend
c     Step 1:
c
c           get solutions of the small matrix.
c
c                btmp holds the initial matrix which is destroyed.
c
c        
         call lsolve(htmp,btmp,resid,resid(1+lwork,1),nend,newrhs,
     1               maxvec,lwork)
         if(prnt(5)) then
            title='solutions of small matrix iteration = '//itoc(iter) 
            call prntrm(title,btmp,nend,newrhs,maxvec,newrhs,iout)
         endif
c      
c     Step 2:
c
c           form the residuals and check for convergence.
c           t contains the all the solution vectors.
c
c
c           when lares is finished, the converged solutions are stored
c           on the disk and the unconverged residuals and solutions appear
c           first in the list.  the b and btmp matrices are updated to
c           reflect this and the current value of newrhs changed to the
c           new number.       
c 
         call lares(pvec,hpvec,b,btmp,rhs,energy,cnverg,resid,
     1              maxerr,t,list,n,nend,newrhs,con,uncon,maxvec,
     2              iter,prnt(5))
         newrhs=uncon
         totrhs=totrhs+con
         if(totrhs.eq.nrhs) then
c
c           all solutions are converged.  copy them in to rhs and quit
c
            return
         else
c
c           all solutions are not converged.  set the error to the largest
c           current error and continue the iteration sequence 
c
            error=min(error,maxerr)
c
c           scale the residuals to get the next set of vectors to
c           be added to the krylov sequence.  leave them in resid
c           for the moment.
c 
            call lavec(energy,diag,resid,n,uncon,iter,prnt(9))
c
c           orthogonalize the scaled residuals to the old vectors
c           to get an additional nout set of linearly independent
c           vectors.
c        
            call abschm(pvec,resid,thresh,n,nend,uncon,nout,
     1                  .true.,.false.)
c
c           orthonormalize these vectors and add them to the
c           existing set.
c
            call gschmt(resid,thresh,n,1,nout,nfinal,.true.,.false.)
c               
c           how many new vectors could be added in principle
c
            numnew = maxvec - nend
c                                             
c           how many will we add
c
            addvec = min(numnew,nfinal)
            write(iout,5) nfinal, addvec
c
c           add them
c
            nbeg = nend + 1            
            nend=nend+addvec
            if(nend.ge.maxvec.and.addvec.eq.0) then
               write(iout,6)
               stop
            endif
            call copy(resid,pvec(1,nbeg),n*addvec)
            if(prnt(10)) then
               call tstovl(pvec,n,nend)
            endif
c
c           we now have either a reconstituted set of vectors or
c           an augmented set of vectors.               
c
c              operate with hamiltonian on these vectors
c
               title='h on initial vectors'
               call honv(ibuf,hbuf,diag,pvec(1,nbeg),hpvec(1,nbeg),
     1                   n,addvec,lenbuf,ntot,incore,iter,hamno,prnt(3))
c
c           update the curent small hamiltonian matrix and right hand sides.
c            
               call bsmall(h,htmp,b,btmp,pvec,hpvec,rhs,energy,newrhs,
     1                     n,nbeg,nend,maxvec,'fill')       
         endif            
      enddo
 99   continue                       
      return
 1    format(/,1x,'davidson linear system solver',/,10x,
     1            'size of matrix                = ',i5,/,10x,
     2            'number of right hand sides    = ',i4,/,10x,
     3            'maximum number of iterations  = ',i4,/,10x,
     4            'maximum number of vectors     = ',i4,/,10x,
     5            'convergence criterion for rms = ',e15.8)      
 2    format(/,5x,'cannot even begin davidson calculation:',/,5x,
     1            'orthonormalization of initial vectors yields null'
     2            '  set')
 3    format(/,5x,'beginning iterations:',/,5x,
     1            'error initialized at ',e15.8) 
 4    format(/,5x,'iteration            = ',i5,
     1       /,5x,'size of vector space = ',i5)
 5    format(/,5x,'number of linearly independent vectors after '
     1            'schmidt step = ',i4,/,5x,'number actually added '
     2            '             = ',i4)
 6    format(/,5x,'size of allowed vector space exceeded. quit')
      end       





