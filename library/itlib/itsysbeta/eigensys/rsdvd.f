*deck rsdvd.f
c***begin prologue     rsdvd
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative eigenvalue
c***author             schneider, barry (nsf)
c***source             
c***purpose            find the eigenvalues of a real symmetric matrix
c***                   having diagonal dominance using the davidson 
c***                   algorithm.  
c***
c***description        the davidson algorithm as programmed here is
c***                   capable of getting the lowest m roots or of
c***                   climbing the eigenvalue tree in a piecewise fashion
c***                   by getting the first m1 roots, then the next m2 roots
c***                   then the next.... roots until all are obtained.
c***                   it can often happen that roots will be missed, 
c***                   especially if one tries to compute too many 
c***                   at each pass.  so, the routine must be used 
c***                   intelligently.  it is typically simple to see if a
c***                   root has been missed by inspection.  the method avoids
c***                   diagonalizing a large matrix for the higher roots by
c***                   explicitly orthogonalizing the guesses to the converged
c***                   eigenvectors of the lower roots.
c***
c***     variable      type       size             description
c***     _______       ____       ____             ___________
c***      ibuf          integer 2*lenbuf      buffer to hold non-zero values of
c***                                          hamiltonian indices
c***      rbuf          real    lenbuf        buffer to hold non-zero values of
c***                                          hamiltonian values
c***      diag          real    n             diagonals of hamiltonian

c***      pvec          real    n*maxvec      input as trials but become
c***                                          the davidson vectors as the
c***                                          calculation proceeds
c***                                          can be same storage as trials
c***      hpvec         real    n*maxvec      action of the hamiltonian on pvec
c***      b, btmp       real    maxvec*maxvec small hamiltonian matrix and a 
c***                                          copy
c***      etmp          real    maxvec        eigenvalues of small matrix
c***      work          real    5*maxvec      scratch array
c***      svec          real    maxvec*maxvec small scratch matrix
c***      resid         real    n*maxvec      residual vectors
c***      eig           real    nroots        converged eigenvalues
c***      rep           real                  nuclear repulsion
c***      cnverg        real                  convergence criterion for a root
c***      thresh        real                  overlap tolerance for accepting 
c***                                          a new davidson vector
c***      n             integer               matrix size
c***      nroots        integer               number of roots to get
c***      ntrials       integer               number of trial vectors
c***      nattim        integer               number of roots to converge at a 
c***                                          pass
c***      maxit         integer               maximum size of davidson space
c***      maxvec        integer               maximum number of vectors 
c***                                          available for storage
c***      lenbuf        integer               size of buffer for non-zero 
c***                                          labels and matrix elements
c***      ntot          integer               total number of non-zero
c***                                          off-diagonal matrix elements
c***                                          of matrix  
c***      header        character             labels for reading buffers
c***      prnt          logical               print flags   
c***references         
c
c***routines called    
c***end prologue       rsdvd
      subroutine rsdvd(ibuf,rbuf,diag,pvec,hpvec,b,btmp,etmp,work,
     1                 svec,resid,eig,rep,cnverg,thresh,n,nroots,
     2                 ntrials,nattim,maxit,maxvec,lenbuf,ntot,
     3                 prnt)
      implicit integer (a-z)
      real*8 rbuf, diag, pvec, hpvec, svec, work
      real*8 b, btmp, etmp, eig
      real*8 resid, error, cnverg, thresh, tmp, rep
      real*8 maxerr, zero, one, nrzero, eci
      logical prnt, incore
      character*8 cntrl
      character*5 itoc
      character*80 title
      dimension rbuf(lenbuf), ibuf(2,lenbuf), diag(n)
      dimension pvec(n,maxvec), hpvec(n,maxvec)
      dimension b(maxvec,maxvec), btmp(maxvec,maxvec)
      dimension etmp(maxvec), svec(maxvec,*)
      dimension eig(nroots)
      dimension resid(n,maxvec), work(5*maxvec)
      dimension prnt(11)
c
c     only common block.  it contains the input and output units.
c     typically inp=5 and iout=6.
c 
      common/io/inp, iout
      data zero, one/ 0.d0, 1.d0/
      data nrzero / 1.0d-06 /
      data incore/.false./
      data file/90/
      open(unit=70,file='bliu',err=99,form='unformatted',
     1     status='unknown')
      open(unit=80,file='davidson',err=99,form='unformatted',
     1     status='unknown')
c
c     if ntot is small enough, the non-zero elements and indices will be
c     assumed to reside in the buffers and not be read from disk over and
c     over again. 
c
      if(ntot.le.lenbuf) then
         incore=.true.
	 write(iout,*) 'Matrix Vector Multiply in Core'
      endif
      write(iout,1) nroots, nattim, ntrials, maxit, maxvec, cnverg
c
c     put the trials on the bliu scratch file
c
      do 2000 i=1,ntrials
         write(70) (pvec(j,i), j=1,n)
 2000 continue   
      if(prnt(1)) then
         title='trials'
         call prntfm(title,pvec,n,ntrials,n,ntrials,iout)
      endif
c
c     find out how many passes are needed to get all the roots
c      
      npass=nroots/nattim
      nleft=nroots-npass*nattim
      if(nleft.ne.0) then
         npass=npass+1
      else
         nleft=nattim
      endif            
      count=0
      root0=0
      rootn=0
      num2do=nattim
      numvec=ntrials
      do 1000 ipass=1,npass
         if(ipass.eq.npass) then
            num2do=nleft
         endif
         root0=rootn+1
         rootn=rootn+num2do    
         write(iout,2) ipass, root0, rootn      
c-----------------------------------------------------------------------c
c                                                                       c
c                    Initialization Section                             c
c                                                                       c
c-----------------------------------------------------------------------c
c
c        get the initial set of vectors by orthogonalizing the trials
c        to the set of roots which have been converged.  this will ensure
c        that the subspace is orthogonal to previously converged vectors. 
c
         num=numvec
         if(count.ne.0) then
            call rdciv(resid,n,count)
            call abschm(resid,pvec,thresh,n,count,num,nout,
     1                  .true.,.false.)
            num=nout
         endif
         if(num.ne.0) then
            call gschmt(pvec,thresh,n,1,num,nout,.true.,.false.)
         endif
         if(nout.eq.0) then
            write(iout,3)
            return
         endif                  
         if(prnt(2)) then
            title='initial vectors'
            call prntfm(title,pvec,n,nout,n,nout,iout)
         endif
c
c        initialize the effect of the hamiltonian on these vectors.
c
         title='h on initial vectors'
         call honv(ibuf,rbuf,diag,pvec(1,1),hpvec(1,1),
     1             n,nout,lenbuf,ntot,incore,0,file,prnt(3))
c
c        initialize the small hamiltonian matrix.
c            
         call hsmall(b,btmp,pvec,hpvec,n,ndum,nout,maxvec,'initialize')
c----------------------------------------------------------------------c
c                                                                      c
c                    Iteration Sequence                                c
c                                                                      c
c     iteration is continued until all of the roots are converged      c 
c     or if convergence is not achieved some preset maximum number of  c
c     iterations are performed.                                        c
c                                                                      c
c----------------------------------------------------------------------c
         nend=nout 
         iter=0
         error=1.d+10
         write(iout,4) error
         cntrl='continue'
         do while ( cntrl.eq.'continue'.and.iter.lt.maxit )
            iter = iter + 1
c        Step 1:
c
c           get eigenvalues and eigenvectors of the small matrix.
c
c                btmp holds the initial matrix which is destroyed.
c                svec has the transformation matrix.
c                note that resid is used as temporary storage in rdiag
c                and vscale.
c
c
            if(prnt(4)) then
               title='initial small matrix iteration = '//itoc(iter)
               call prntfm(title,btmp,nend,nend,maxvec,maxvec,iout)
            endif
c        
            call rdiag(btmp,etmp,svec,work,rep,iter,maxvec,
     1                 nend,prnt(5))
c      
c        Step 2:
c
c           transform pvec and hpvec to the new basis defined by the
c           diagonalization in rdiag.  the small matrix becomes diagonal under
c           this transformation and we fill it with the eigenvalues.
c           then form the residuals and check for convergence.  the converged
c           eigenvalues are stored and the eigenpairs written to disk.  the
c           unconverged residuals are moved from their current positions in resid
c           to the leading positions.
c
            write(iout,5) iter, nend
            call frmres(pvec,hpvec,svec,eig,etmp,b,btmp,rep,cnverg,
     1                  resid,maxerr,n,nend,num2do,con,uncon,
     2                  maxvec,iter,prnt(5),count)
c
c     
c           check to see if all num2do roots are converged
c           if so, we are done and can print results and then proceed to 
c           the next set of roots.
c
            if(con.eq.num2do) then
c
               nlft=min(nend-num2do,ntrials)
               remain=ntrials-nlft
               write(iout,6)
               do 10 i=1,num2do
                  actual=i+count
                  write(iout,7) actual, eig(actual)+rep
 10            continue
c
c              begin next set of roots with the best available
c              vectors.  these are the the remaining vectors coming 
c              from the diagonalization of the small matrix augmented
c              with the remaining trials.
               if(actual.lt.nroots) then
                  write(iout,8) nlft, remain
                  if(nlft.gt.0) then
                     write(iout,*) ' Doing the vectors moves'
                     call copy(pvec(1,num2do+1),pvec,n*nlft)  
                     if(remain.gt.0) then
c                       
c                       skip over the trials not needed
c
                        rewind(70)
                        do 2001 skp=1,nlft
                           read(70) (pvec(jj,nlft+1), jj=1,n)
 2001                   continue
c
c                       read in the trials needed
c   
                        cnt=nlft
                        do 2002 addin=1, remain   
                           cnt=cnt+1
                           read(70) (pvec(jj,cnt),jj=1,n)
 2002                   continue   
                     endif
                  else
                     rewind(70)
                     do 2003 addin=1,ntrials
                         read(70) (pvec(jj,addin), jj=1,n)
 2003                continue
                  endif
               endif
               numvec=ntrials
               cntrl='finished'
            else
c
c           all roots are not converged.  set the error to the largest
c           current error and either restart or continue the 
c           iteration sequence.
c
               error=min(error,maxerr)
c
c           how many new vectors could be added in principle
c
               numnew = maxvec - nend
c                                             
c           how many will we add
c
               addvec = min(numnew,uncon)
c              
c           check if the number of old plus new vectors will exceed
c           the maxvec limit to determine if a restart is needed.
c
               chkno = nend + addvec    
               if(chkno.ge.maxvec) then
c             
c              the maximum number of vectors is exceeded.
c
                  nend = min(2*num2do,nend)
                  write(iout,9) nend
                  numnew = maxvec - nend
                  addvec = min(numnew,uncon)
               endif
               write(iout,11) addvec
c
c              maximum number of vectors is still within the allowed
c              limits.  add vectors to the set from the unconverged
c              residuals and put them after the current vectors.
c
               if(addvec.ne.0) then
                  oldnum = nend   
                  nbeg = nend + 1
                  call newvec(etmp,diag,resid,pvec(1,nbeg),
     1                        n,addvec,iter,prnt(9))
c
c              orthonormalize the new trials to the converged roots
c              and then to the old vectors to get an additional 
c              nout vectors.
c        
                  num=addvec
                  nout=num
                  if(count.ne.0) then
                     call rdciv(resid,n,count)
                     call abschm(resid,pvec(1,nbeg),thresh,n,count,
     1                           addvec,num,.true.,.false.)
                  endif
                  nend=nend+num
                  if(num.ne.0) then
                     call gschmt(pvec,thresh,n,nbeg,nend,
     1                           nout,.true.,.false.)
                  endif
                  if(num.eq.0.or.nout.eq.0) then
                     write(iout,12) num, nout
                    return
                  endif
c                
                  if(nout.eq.0) then
c
c                 no more vectors write out unconverged results.
c               
                     write(iout,13)
                     return
                  endif   
                  nend=nbeg+nout-1
                  if(prnt(10)) then
                     call tstovl(pvec,n,nend)
                  endif
               endif            
c
c           we now have either a reconstituted set of vectors or
c           an augmented set of vectors.               
c
c              operate with hamiltonian on these vectors
c
               call honv(ibuf,rbuf,diag,pvec(1,nbeg),hpvec(1,nbeg),
     1                   n,nout,lenbuf,ntot,incore,iter,file,prnt(3))
c
c           update the curent small hamiltonian matrix.
c            
               call hsmall(b,btmp,pvec,hpvec,n,nbeg,nend,maxvec,'fill')
            endif
         enddo
         if(iter.ge.maxit) then
            write(iout,14)   
            return
         endif        
         count=count+num2do
 1000 continue
      write(iout,15)
      do 20 i=1,nroots
         eci=eig(i)+rep
         write(iout,7) i, eci
 20   continue           
      return
 99   write(iout,*) 'error in opening unit 70'
 1    format(/,1x,'davidson eigenvalue solver',/,10x,
     1            'number of roots               = ',i4,/,10x,
     2            'number of roots at a time     = ',i4,/,10x,
     3            'number of trials              = ',i4,/,10x,
     4            'maximum number of iterations  = ',i4,/,10x,
     5            'maximum number of vectors     = ',i4,/,10x,
     6            'convergence criterion for rms = ',e15.8)      
 2    format(/,5x,'pass = ',i3,/,5x,
     1            'processing root = ',i4,' to root =',i4)
 3    format(/,5x,'cannot even begin davidson calculation:',/,5x,
     1            'orthonormalization of initial vectors yields null '
     2            'set')      
 4    format(/,5x,'beginning davidson iterations:',/,5x,
     1            'initial error = 'e15.8) 
 5    format (/,5x,'cycle = ',i4,2x,'no. vectors = ',i4)
 6    format(/,30x,'summary',/,17x,'root',22x,'energy')
 7    format(15x,i4,20x,f15.8)
 8    format(/,5x,'number of added vectors from unconverged roots = ',
     1                                                        i5,/,5x,
     2            'number of added vectors from original trial set = ',
     3                                                        i5)
 9    format(/,10x,'***** maximum number of vectors exceeded *****',
     1        /10x,'      contract back to ',1x,i4,' vectors')
 11   format(/,5x,'number of vectors added = ',i4)                   
 12   format(/,5x,'cannot continue',/,5x,
     1            'num = ',i5,2x,'nout = ',i5,/,5x,
     2            'write out results and return')
 13   format(/1x,'no more orthonormal vectors can be added.',/,5x,
     1           'quit and return unconverged results')
 14   format(/,5x,'iteration limit exceeded.',/,5x,
     1            'quit and return to main')
 15   format(///,30x,'final summary',/,17x,'root',22x,'energy')

      end       










