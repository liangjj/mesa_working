*deck rpa.f
c***begin prologue     rpa
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative eigenvalue rpa
c***author             schneider, barry (nsf)
c***source             
c***purpose            iterative davidson routine for rpa 
c***                   eigenvalue problem.
c***                   
c***references         
c
c***routines called    
c***end prologue       rpa
      subroutine rpa(hbufa,ibufa,diaga,hbufb,ibufb,trials,etrial,
     1               vecl,vecr,psi,pvec,hpvec,vec,b,btmp,etmp,
     2               resid,er,ei,t1,t2,lwork,cnverg,thresh,
     3               n,nroots,ntrials,nattim,maxit,maxvec,
     4               lenbuf,drctv,incore,neapb,neamb,prnt,mattyp)
      implicit integer (a-z)
      real*8 hbufa, diaga, hbufb, trials, etrial, vecl, vecr, psi
      real*8 pvec, hpvec, vec, b, btmp, etmp, resid, er, ei
      real*8 t1, t2, error, cnverg, thresh, sdot, tmp, maxerr
      real*8 scale
      logical prnt, drctv, incore
      character*(*) mattyp
      character*5 itoc
      character*80 title
      dimension hbufa(lenbuf), hbufb(lenbuf), ibufa(2,lenbuf)
      dimension ibufb(2,lenbuf), diaga(n), trials(n,ntrials)
      dimension etrial(ntrials), vecl(maxvec,*), vecr(maxvec,*)
      dimension pvec(n,maxvec), hpvec(n,maxvec), vec(n,maxvec)
      dimension b(maxvec,maxvec), btmp(maxvec,maxvec), etmp(maxvec)
      dimension resid(n,maxvec), er(maxvec), psi(n,*)
      dimension ei(maxvec), t1(n,*), t2(n,*)
      dimension prnt(11)
      common/io/inp, iout
      scale=1.d0
      if(.not.incore) then
         call iosys('read integer "number of a plus b matrix '//
     1              'elements" from ham',1,neapb,0,' ')
         call iosys('read integer "number of a minus b matrix '//
     1              'elements" from ham',1,neamb,0,' ')
      endif
      if(drctv) then
         write(iout,1) nroots, ntrials, maxit, maxvec, cnverg
      else
         write(iout,2) nroots, ntrials, maxit, maxvec, cnverg
      endif
      if(prnt(1)) then
         title='trials'
         call prntfm(title,trials,n,ntrials,n,ntrials,iout)
      endif
c-----------------------------------------------------------------------c
c                                                                       c
c                    Initialization Section                             c
c                                                                       c
c-----------------------------------------------------------------------c
c
c     initialize the set of input vectors as the 
c     orthonormalized set of trial vectors
c
      call copy(trials,pvec,n*ntrials)
      nbeg=1
      nend=ntrials
      call gschmt(pvec,thresh,n,1,nend,nout,.true.,.false.)
      if(nend.eq.0) then
         write(iout,4)
         call lnkerr('quit rpa davidson')
      else
         nend=nout
         if(prnt(2)) then
            title='initial vectors'
            call prntfm(title,pvec(1,nbeg),n,nout,n,nout,iout)
         endif
      endif
c
c    initialize the effect of the hamiltonian on these vectors.
c
      title='h on initial vectors'
      call rpaonv(hbufa,hbufb,ibufa,ibufb,diaga,pvec(1,nbeg),
     1            hpvec(1,nbeg),t1,n,nout,lenbuf,neapb,neamb,
     2            incore,title,0,prnt(3))
c
c     initialize the small hamiltonian matrix.
c            
      do 10 i=1,nend
         do 20 j=1,i
            b(i,j) = sdot(n,pvec(1,i),1,hpvec(1,j),1)
            btmp(i,j) = b(i,j)
            b(j,i) = sdot(n,pvec(1,j),1,hpvec(1,i),1)
            btmp(j,i) = b(j,i)
 20      continue   
 10   continue
c----------------------------------------------------------------------c
c                                                                      c
c                    Iteration Sequence                                c
c                                                                      c
c     iteration is continued until all of the roots are converged      c 
c     or if convergence is not achieved some preset maximum number of  c
c     iterations are performed.                                        c
c                                                                      c
c----------------------------------------------------------------------c
      iter=0
      error=1.d+10
      write(iout,5) error
      do while ( error.gt.cnverg.and.iter.lt.maxit )
         iter = iter + 1
         write(iout,6) iter, nend
c     Step 1:
c
c           get eigenvalues and eigenvectors of the small matrix.
c
c                btmp holds the initial matrix which is destroyed.
c                vec has the transformation matrix.
c                note that resid is used as temporary storage in rdiag
c                and vscale.
c
         title='iteration = '//itoc(iter)//' diagonalizing  matrix '
     1                      //'of dimension N = '//itoc(nend)
         write(iout,7) title
c        
         if(prnt(4)) then
            title='small matrix'
            call prntfm(title,btmp,nend,nend,maxvec,maxvec,iout)
         endif
c         call rpadiag(btmp,er,ei,vecl,vecr,t1,maxvec,nend,lwork)
         call rpadiag(btmp,etmp,ei,vecl,vecr,t1,maxvec,nend,lwork)
         if(prnt(5)) then
            title='eigenvalues of small matrix iteration = '//itoc(iter) 
c            call prntfm(title,er,nend,1,nend,1,iout)
            call prntfm(title,etmp,nend,1,nend,1,iout)
         endif
c      
c     Step 2:
c
c           form the residuals and check for convergence.
c           t1 contains the transformed vectors and t2 the transformed
c           hamiltonian on vectors.
c
c
c        unconverged residuals are in resid and the corresponding eigenvalues
c        are in etmp.
c 
c         call copy(er,etmp,nroots)              
         call frmres(psi,pvec,hpvec,vecr,er,etmp,scale,cnverg,resid,
     1               maxerr,t1,t2,n,nend,nroots,con,
     2               uncon,maxvec,iter,prnt(5))
         if(con.eq.nroots) then
c
c           all roots are converged.  lets quit
            write(iout,8)
            do 30 i=1,nroots
               write(iout,9) i, er(i)
 30         continue
            return
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
c            if(addvec.eq.0) then
c               write(iout,13)
c               return
c            endif
c              
c           check if the number of old plus new vectors will exceed
c           the maxvec limit to determine if a restart is needed.
c
            chkno = nend + addvec    
            if(chkno.ge.maxvec) then
c             
c              the maximum number of vectors is exceeded.
c
c              transform the current set of nend vectors to a new
c              set of orthonormal vectors based on the right
c              eigenvectors of the current small matrix.  that requires
c              a schmidt process applied to the right eigenvector
c              matrix.  once that is done transform the operations of
c              the hamiltonian on the old basis to that on the new
c              basis and transform the small matrix to the new basis.
c              the small matrix will be upper triangular in the new
c              basis with the eigenvalues along the diagonal. 
c
               call rpainit(pvec,hpvec,t1,t2,vecr,vecl,btmp,btmp,b,
     1                      thresh,n,nend,maxvec,prnt(11))
c
c              now take a subset of the transformed vectors plus the
c              added vectors as a new set and re-initiate the process.
c
c               nend = maxvec - addvec
               nend = min(nend,nroots)
               numnew = maxvec - nend
               addvec = min(numnew,uncon)
c               addvec = 0
c               nbeg = 1
               write(iout,11) addvec
c               nout = nend
            endif
c
c              maximum number of vectors is still within the allowed
c              limits.  add vectors to the set from the unconverged
c              residuals and put them after the current vectors.
c
            if(addvec.ne.0) then
               nbeg = nend + 1
               call rpavec(etmp,diaga,resid,pvec(1,nbeg),n,addvec,
c               call rpavec(etmp,etrial,resid,pvec(1,nbeg),n,addvec,
     1                     iter,prnt(9))
               nend=nend+addvec
c
c              orthonormalize the new trials to the old vectors
c              to get an additional nout vectors.
c        
               call gschmt(pvec,thresh,n,nbeg,nend,nout,.true.,.false.)
c               
               if(nout.eq.0) then
c
c                 no more vectors write out unconverged results.
c               
                  write(iout,12)
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
            call rpaonv(hbufa,hbufb,ibufa,ibufb,diaga,pvec(1,nbeg),
     1                  hpvec(1,nbeg),t1,n,nout,lenbuf,neapb,neamb,
     2                  incore,title,iter,prnt(3))
c
c           update the curent small hamiltonian matrix.
c            
            do 40 i=1,nend
               do 50 j=nbeg,nend
                   b(i,j) = sdot(n,pvec(1,i),1,hpvec(1,j),1)
                   b(j,i) = sdot(n,pvec(1,j),1,hpvec(1,i),1)
 50            continue   
 40         continue
            do 60 i=1,nend   
               do 70 j=1,i
                  btmp(i,j) = b(i,j)
                  btmp(j,i) = b(j,i)
 70            continue
 60         continue  
         endif
      enddo        
      call iosys('write integer "size of rpa vector '//
     1           'space" to ham',1,nend,0,' ')
      call iosys('write real "rpa eigenvalues" to ham',
     1            nend,eig,0,' ')
      call iosys('write real "rpa vectors" to ham',
     1            nend*n,pvec,0,' ')            
      return
 1    format(/,1x,'rpa davidson eigenvalue solver using '
     1            'preconditioning',/,10x,
     2            'number of roots               = ',i4,/,10x,
     3            'number of trials              = ',i4,/,10x,
     4            'maximum number of iterations  = ',i4,/,10x,
     5            'maximum number of vectors     = ',i4,/,10x,     
     6            'convergence criterion for rms = ',e15.8)      
 2    format(/,1x,'rpa davidson eigenvalue solver',/,10x,
     1            'number of roots               = ',i4,/,10x,
     2            'number of trials              = ',i4,/,10x,
     3            'maximum number of iterations  = ',i4,/,10x,
     4            'maximum number of vectors     = ',i4,/,10x,
     5            'convergence criterion for rms = ',e15.8)      
 4    format(/,5x,'cannot even begin davidson calculation:',/,5x,
     1            'orthonormalization of initial vectors yields null'
     2            '  set')
 5    format(/,5x,'beginning davidson iterations:',/,5x,
     1            'error initialized at ',e15.8) 
 6    format(/,1x,'beginning next davidson cycle',/,5x,
     1             'iteration = ',i4,1x,'size of vector space = ',i5)
 7    format(/,5x,a80)
 8    format(/,30x,'summary',/,17x,'root',22x,'energy')
 9    format(15x,i4,20x,e15.8)
 11   format(/,5x,'number of added vectors will exceed maxvec',/,5x,
     1            'number of vectors actually added = ',i4)                   
 12   format(/1x,'no more orthonormal vectors can be added')
 13   format(/,1x,'no more vectors can be added')
      end       


