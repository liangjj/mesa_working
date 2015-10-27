*deck rsdvd.f
c***begin prologue     rsdvd
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative eigenvalue
c***author             schneider, barry (nsf)
c***source             
c***purpose            iterative davidson eigenvalue.
c***                   
c***references         
c
c***routines called    
c***end prologue       rsdvd
      subroutine rsdvd(hbuf,ibuf,diag,trials,eig01,eig02,eig03,
     1                 u01,u02,u03,psi,pvec,hpvec,vec,b,btmp,etmp,
     2                 svec,resid,eig,t1,t2,scale,cnverg,thresh,
     3                 n1,n2,n3,dim,n,nroots,ntrials,nattim,
     4                 maxit,maxvec,lenbuf,drctv,incore,ntot,
     5                 prnt,mattyp)
      implicit integer (a-z)
      real*8 hbuf, diag, trials, pvec, hpvec, psi, vec, svec
      real*8 b, btmp, etmp, eig, eig01, eig02, eig03, u01, u02, u03
      real*8 resid, t1, t2, error, cnverg, thresh, sdot, scale, tmp
      real*8 maxerr, zero, one, nrzero
      complex*16 cdum
      logical prnt, drctv, incore
      character*(*) mattyp
      character*5 itoc
      character*8 cntrl
      character*1 it
      character*80 title
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n), trials(n,ntrials)
      dimension pvec(n,maxvec), hpvec(n,maxvec), psi(n,nroots)
      dimension vec(n,maxvec), b(maxvec,maxvec), btmp(maxvec,maxvec)
      dimension etmp(maxvec), svec(maxvec,*)
      dimension eig(nroots), eig01(n1), eig02(n2), eig03(n3)
      dimension u01(n1,n1), u02(n2,n2), u03(n3,n3)
      dimension resid(n,maxvec), t1(n,*), t2(n,*)
      dimension prnt(11)
      real*8 ham, temp
      dimension ham(100,100), temp(100,100), ipvt(100)
      common/io/inp, iout
      data zero, one/ 0.d0, 1.d0/
      data nrzero / 1.0d-06 /
      it=itoc(dim)
      if(.not.incore) then
         call iosys('read integer "'//it//
     1              'd number of elements" from ham',1,ntot,0,' ')
      endif
      if(drctv) then
         write(iout,1) nroots, nattim, ntrials, maxit, maxvec, cnverg
      else
         write(iout,2) nroots, nattim, ntrials, maxit, maxvec, cnverg
      endif
      if(prnt(1)) then
         title='trials'
         call prntfm(title,trials,n,ntrials,n,ntrials,iout)
      endif
      call copy(trials,pvec,ntrials*n)
      write(iout,3) scale
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
         write(iout,4) ipass, root0, rootn      
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
            call abschm(psi,pvec,thresh,n,count,num,nout,.true.)
            num=nout
         endif
         if(num.ne.0) then
            call gschmt(pvec,thresh,n,1,num,nout,.true.)
         endif
         if(nout.eq.0) then
            write(iout,5)
            call lnkerr('quit davidson. no more trial vectors '//
     1                  'possible')
         endif                  
         if(prnt(2)) then
            title='initial vectors'
            call prntfm(title,pvec,n,nout,n,nout,iout)
         endif
c
c        initialize the effect of the hamiltonian on these vectors.
c
         title='h on initial vectors'
         call honv(hbuf,hbuf,ibuf,diag,diag,pvec(1,1),
     1             pvec(1,1),hpvec(1,1),hpvec(1,1),
     2             n,nout,lenbuf,ntot,it,incore,mattyp,
     3             title,0,prnt(3))
c
c        initialize the small hamiltonian matrix.
c            
         do 10 i=1,nout
            do 20 j=1,i
               b(i,j) = sdot(n,pvec(1,i),1,hpvec(1,j),1)
               btmp(i,j) = b(i,j)
               b(j,i) = b(i,j)
               btmp(j,i) = b(j,i)
 20         continue   
 10      continue
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
         write(iout,6) error
         cntrl='continue'
         do while ( cntrl.eq.'continue'.and.iter.lt.maxit )
            iter = iter + 1
            write(iout,7)
c        Step 1:
c
c           get eigenvalues and eigenvectors of the small matrix.
c
c                btmp holds the initial matrix which is destroyed.
c                svec has the transformation matrix.
c                note that resid is used as temporary storage in rdiag
c                and vscale.
c
            title='iteration = '//itoc(iter)//' diagonalizing '//
     1            'matrix of dimension N = '//itoc(nend)
            write(iout,8) title
            if(prnt(4)) then
               title='initial small matrix'
               call prntfm(title,btmp,nend,nend,maxvec,maxvec,iout)
            endif
c        
            call rdiag(btmp,btmp,etmp,etmp,svec,svec,svec,
     1                 resid,ldum,resid,maxvec,nend,mattyp)
            call vscale(resid,etmp,scale,nend)
            if(prnt(5)) then
               title='eigenvalues of small matrix iteration = '
     1               //itoc(iter)
               call prntfm(title,resid,nend,1,nend,1,iout)
            endif
c      
c        Step 2:
c
c           form the residuals and check for convergence.
c           t1 contains the transformed vectors and t2 the transformed
c           hamiltonian on vectors.
c
c
c          unconverged residuals are in resid and the corresponding 
c          eigenvalues are in etmp.
c 
            call frmres(psi(1,count+1),pvec,hpvec,svec,eig(count+1),
     1                  etmp,scale,cnverg,resid,maxerr,t1,t2,n,nend,
     2                  num2do,con,uncon,maxvec,iter,prnt(5))
c
c           note that t1 contains the best set of vectors at this point
c           in the calculation.
c     
            if(con.eq.num2do) then
c
c           all num2do roots are converged.  print results and then 
c           proceed to next set of roots
               nlft=min(nend-num2do,ntrials)
               remain=ntrials-nlft
               write(iout,9)
               do 30 i=1,num2do
                  actual=i+count
                  write(iout,11) actual, eig(actual)
 30            continue
c
c              begin next set of roots with the best available
c              vectors.  these are the converged roots plus the remaining
c              vectors coming from the diagonalization of the small matrix
c              if there are not enough vectors, fill out the remaining
c              subspace with the trials.
               if(actual.lt.nroots) then
                  write(iout,12) nlft, remain
                  if(nlft.gt.0) then
                     call copy(t1(1,num2do+1),pvec,n*nlft)  
                     if(remain.gt.0) then
                        call copy(trials(1,nlft+1),pvec(1,nlft+1),
     1                            n*remain)
                     endif
                  else
                     call copy(trials,pvec,n*ntrials)
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
c              put the current best vectors and the operation of the
c              hamiltonian on these vectors in the proper positions,
c              transform the small matrix to the new basis
c              and reset the counters.
c
                  call reinit(t1,t2,pvec,hpvec,svec,b,btmp,n,nend,
     1                        num2do,maxvec,prnt(11))
                  nend = min(2*num2do,nend)
                  write(iout,18) nend
                  numnew = maxvec - nend
                  addvec = min(numnew,uncon)
               endif
               write(iout,13) addvec
c
c              maximum number of vectors is still within the allowed
c              limits.  add vectors to the set from the unconverged
c              residuals and put them after the current vectors.
c
               if(addvec.ne.0) then
                  oldnum = nend   
                  nbeg = nend + 1
                  call newvec(etmp,diag,resid,pvec(1,nbeg),
     1                        eig01,eig02,eig03,u01,u02,u03,
     2                        t1,t2,dim,n1,n2,n3,n,addvec,drctv,
     3                        iter,prnt(9))
c                  call exvec(ham,etmp,hbuf,ibuf,diag,resid,
c     1                       pvec(1,nbeg),temp,ipvt,n,addvec,
c     2                       ntot,iter,.true.)
c
c              orthonormalize the new trials to the converged roots
c              and then to the old vectors to get an additional 
c              nout vectors.
c        
                  num=addvec
                  nout=num
                  if(count.ne.0) then
                     call abschm(psi,pvec(1,nbeg),thresh,n,count,
     1                           addvec,num,.true.)
                  endif
                  nend=nend+num
                  if(num.ne.0) then
                     call gschmt(pvec,thresh,n,nbeg,nend,
     1                           nout,.true.)
                  endif
                  if(num.eq.0.or.nout.eq.0) then
                     write(iout,14)
                     call iosys('write integer "size of '//
     1                          'davidson vector space" to ham',
     2                           1,nend,0,' ')
                    call vscale(eig,eig,scale,nend)
                    call iosys('write real "davidson eigenvalues" '//
     1                         'to ham', nend,eig,0,' ')
                    call iosys('write real "davidson vectors" '//
     1                         'to ham', nend*n,pvec,0,' ')            
                    return
                  endif
c                
                  if(nout.eq.0) then
c
c                 no more vectors write out unconverged results.
c               
                     write(iout,15)
                     return
                  endif   
                  nend=nbeg+nout-1
                  if(prnt(10)) then
                     call tstovl(pvec,pvec,n,nend,mattyp)
                  endif
               endif            
c
c           we now have either a reconstituted set of vectors or
c           an augmented set of vectors.               
c
c              operate with hamiltonian on these vectors
c
               title='h on initial vectors'
               call honv(hbuf,hbuf,ibuf,diag,diag,pvec(1,nbeg),
     1                   pvec(1,nbeg),hpvec(1,nbeg),hpvec(1,nbeg),
     2                   n,nout,lenbuf,ntot,it,incore,mattyp,
     3                   title,iter,prnt(3))
c
c           update the curent small hamiltonian matrix.
c            
               do 40 i=1,nend
                  do 50 j=nbeg,nend
                     b(i,j) = sdot(n,pvec(1,i),1,hpvec(1,j),1)
                     b(j,i) = b(i,j)
 50               continue   
 40            continue
               do 60 i=1,nend   
                  do 70 j=1,i
                     btmp(i,j) = b(i,j)
                     btmp(j,i) = b(i,j)
 70               continue
 60            continue  
            endif
         enddo
         if(iter.ge.maxit) then
            write(iout,16)   
            return
         endif        
         count=count+num2do
 1000 continue
      write(iout,17)
      do 80 i=1,nroots
         write(iout,11) i, eig(i)
 80   continue           
      call iosys('write integer "size of davidson vector '//
     1           'space" to ham',1,nend,0,' ')
      call vscale(eig,eig,scale,nend)
      call iosys('write real "davidson eigenvalues" to ham',
     1            nend,eig,0,' ')
      call iosys('write real "davidson vectors" to ham',
     1            nend*n,pvec,0,' ')            
      return
 1    format(/,1x,'davidson eigenvalue solver using preconditioning',
     1                                                  /,10x,
     2            'number of roots               = ',i4,/,10x,
     3            'number of roots at a time     = ',i4,/,10x,
     4            'number of trials              = ',i4,/,10x,     
     5            'maximum number of iterations  = ',i4,/,10x,
     6            'maximum number of vectors     = ',i4,/,10x,     
     7            'convergence criterion for rms = ',e15.8)      
 2    format(/,1x,'davidson eigenvalue solver',/,10x,
     1            'number of roots               = ',i4,/,10x,
     2            'number of roots at a time     = ',i4,/,10x,
     3            'number of trials              = ',i4,/,10x,
     4            'maximum number of iterations  = ',i4,/,10x,
     5            'maximum number of vectors     = ',i4,/,10x,
     6            'convergence criterion for rms = ',e15.8)      
 3    format(/,5x,'hamiltonian scale factor = ',e15.8)
 4    format(/,5x,'pass = ',i3,/,5x,
     1            'processing root = ',i4,' to root =',i4)
 5    format(/,5x,'cannot even begin davidson calculation:',/,5x,
     1            'orthonormalization of initial vectors yields null '
     2            'set')      
 6    format(/,5x,'beginning davidson iterations:',/,5x,
     1            'error initialized at ',e15.8) 
 7    format(/,5x,'beginning next davidson cycle')
 8    format(/,5x,a80)
 9    format(/,30x,'summary',/,17x,'root',22x,'energy')
 11   format(15x,i4,20x,e15.8)
 12   format(/,5x,'number of added vectors from unconverged roots = ',
     1                                                        i5,/,5x,
     2            'number of added vectors from original trial set = ',
     3                                                        i5)
 13   format(/,5x,'number of added vectors will exceed maxvec',/,5x,
     1            'number of vectors actually added = ',i4)                   
 14   format(/,5x,'cannot continue.  write out results and return')
 15   format(/1x,'no more orthonormal vectors can be added')
 16   format(/,5x,'iteration limit exceeded. quit')
 17   format(///,30x,'final summary',/,17x,'root',22x,'energy')
 18   format(/,10x,'***** maximum number of vectors exceeded *****',
     1        /10x,'      contract back to ',1x,i4,' vectors')
      end       


