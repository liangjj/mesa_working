*deck cdvd.f
c***begin prologue     cdvd
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative eigenvalue
c***author             schneider, barry (nsf)
c***source             
c***purpose            iterative davidson eigenvalue for unsymmetric
c***                   real eigenvalue equation.
c***references         
c
c***routines called    
c***end prologue       cdvd
      subroutine cdvd(hbuf,ibuf,diag,trials,eig01,eig02,eig03,
     1                u0l1,u0l2,u0l3,u0r1,u0r2,u0r3,pvec,hpvec,vec,
     2                b,btmp,etmp,svecl,svecr,resid,eig,t1,t2,lwork,
     3                scale,cnverg,thresh,n1,n2,n3,dim,n,nroots,ntrials,
     4                maxit,maxvec,lenbuf,drctv,incore,ntot,prnt,mattyp)
      implicit integer (a-z)
      complex*16 hbuf, diag
      complex*16 trials, pvec, hpvec, vec, svecl, svecr 
      complex*16 b, btmp, etmp, eig, eig01, eig02, eig03
      complex*16 u0l1, u0l2, u0l3, u0r1, u0r2, u0r3
      complex*16 resid, t1, t2, cdotc, tmp
      real*8 maxerr, error, cnverg, thresh, scale
      real*8 zero, one, nrzero
      logical prnt, drctv, incore
      character*(*) mattyp
      character*3 itoc
      character*1 it
      character*80 title
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n), trials(n,ntrials)
      dimension pvec(n,maxvec), hpvec(n,maxvec)
      dimension vec(n,maxvec), b(maxvec,maxvec), btmp(maxvec,maxvec)
      dimension etmp(maxvec), eig(ntrials)
      dimension eig01(n1), eig02(n2), eig03(n3)
      dimension u0l1(n1,n1), u0l2(n2,n2), u0l3(n3,n3)
      dimension u0r1(n1,n1), u0r2(n2,n2), u0r3(n3,n3)
      dimension resid(n,maxvec), t1(n,*), t2(n,*), svecl(maxvec,maxvec)
      dimension svecr(maxvec,maxvec)
      dimension prnt(11)
      common/io/inp, iout
      data zero, one/ 0.d0, 1.d0/
      data nrzero / 1.0d-06 /
      it=itoc(dim)
      if(.not.incore) then
         call iosys('read integer "'//it//
     1              'd number of elements" from ham',1,ntot,0,' ')
      endif
      if(drctv) then
         write(iout,1) nroots, ntrials, maxit, maxvec, cnverg
      else
         write(iout,2) nroots, ntrials, maxit, maxvec, cnverg
      endif
      if(prnt(1)) then
         title='trials'
         call prntcm(title,trials,n,ntrials,n,ntrials,iout)
      endif
      write(iout,3) scale
c-----------------------------------------------------------------------c
c                                                                       c
c                    Initialization Section                             c
c                                                                       c
c-----------------------------------------------------------------------c
c
c     initialize the set of input vectors as the 
c     orthonormalized set of trial vectors
c
      call cc2opy(trials,pvec,n*ntrials)
      nbeg=1
      nend=ntrials
      call cschmt(pvec,thresh,n,1,nend,nout,.true.)
      if(nend.eq.0) then
         write(iout,4)
         call lnkerr('quit davidson')
      else
         nend=nout
         if(prnt(2)) then
            title='initial vectors'
            call prntcm(title,pvec(1,nbeg),n,nout,n,nout,iout)
         endif
      endif
c
c    initialize the effect of the hamiltonian on these vectors.
c
      title='h on initial vectors'
      call honv(hbuf,hbuf,ibuf,diag,diag,pvec(1,nbeg),
     1          pvec(1,nbeg),hpvec(1,nbeg),hpvec(1,nbeg),
     2          n,nout,lenbuf,ntot,it,incore,mattyp,
     3          title,0,prnt(3))
c
c     initialize the small hamiltonian matrix.
c            
      do 10 i=1,nend
         do 20 j=1,i
            b(i,j) = cdotc(n,pvec(1,i),1,hpvec(1,j),1)
            btmp(i,j) = b(i,j)
            b(j,i) = cdotc(n,pvec(1,j),1,hpvec(1,i),1)
            btmp(j,i) = b(j,i)
 20      continue   
 10   continue
      if(prnt(4)) then
         title='initial small matrix'
         call prntcm(title,btmp,nend,nend,maxvec,maxvec,iout)
      endif
c
c----------------------------------------------------------------------c
c                                                                      c
c             Complex Conjugate the Left Eigenvectors to make          c
c             subsequent manipulations easier.                         c
c                                                                      c
c----------------------------------------------------------------------c
      call cvconj(u0l1,u0l1,n1*n1)
      if(dim.gt.1) then
         call cvconj(u0l2,u0l2,n2*n2)
      endif
      if(dim.gt.2) then
         call cvconj(u0l3,u0l3,n3*n3)
      endif
c            
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
c                svec has the transformation matrix.
c                note that resid is used as temporary storage in rdiag
c                and vscale.
c
         title='iteration = '//itoc(iter)//' diagonalizing  matrix '
     1                      //'of dimension N = '//itoc(nend)
         write(iout,7) title
         if(prnt(4)) then
            title='small matrix'
            call prntcm(title,btmp,nend,nend,maxvec,maxvec,iout)
         endif
c
         call rdiag(btmp,btmp,eig,eig,svecl,svecl,svecr,
     1              t1,lwork,resid,maxvec,nend,mattyp)
         call cvscal(resid,eig,scale,nend)
         if(prnt(5)) then
            title='eigenvalues of small matrix iteration = '//itoc(iter) 
            call prntcm(title,resid,nend,1,nend,1,iout)
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
         call crmres(pvec,hpvec,svecr,eig,etmp,scale,cnverg,resid,
     1               maxerr,t1,t2,n,nend,nroots,con,
     2               uncon,maxvec,iter,prnt(5))
         if(con.eq.nroots) then
c
c           all roots are converged.  lets quit
            write(iout,8)
            do 30 i=1,nroots
               write(iout,9) i, eig(i)
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
c              the right eigenvector matrix is in svecr and b contains
c              the current small matrix.
c              btmp and svecl are free for use.
c
               call creinit(pvec,hpvec,t1,t2,svecr,svecl,
     1                      btmp,btmp,b,thresh,n,nend,nroots,
     2                      maxvec,prnt(11))
               nend = min(nroots,nend)
               numnew = maxvec - nend
               addvec = min(numnew,uncon)
            endif
            write(iout,11) addvec
            if(addvec.ne.0) then

c
c              maximum number of vectors is still within the allowed
c              limits.  add vectors to the set from the unconverged
c              residuals and put them after the current vectors.
c
               nbeg = nend + 1
               call cnvec(etmp,diag,resid,pvec(1,nbeg),eig01,eig02,
     1                    eig03,u0l1,u0l2,u0l3,u0r1,u0r2,u0r3,
     2                    t1,t2,dim,n1,n2,n3,n,addvec,drctv,
     3                    iter,prnt(9))
               nend=nend+addvec
c               title='all vectors before schmidt at this iteration'
c               call prntcm(title,pvec,n,nend,maxvec,maxvec,iout) 
c
c              orthonormalize the new trials to the old vectors
c              to get an additional nout vectors.
c        
               call cschmt(pvec,thresh,n,nbeg,nend,nout,.true.)
c               title='all vectors after schmidt at this iteration'
c               call prntcm(title,pvec,n,nend,maxvec,maxvec,iout) 
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
            call honv(hbuf,hbuf,ibuf,diag,diag,pvec(1,nbeg),
     1                pvec(1,nbeg),hpvec(1,nbeg),hpvec(1,nbeg),
     2                n,nout,lenbuf,ntot,it,incore,mattyp,
     3                title,iter,prnt(3))
c
c           update the curent small hamiltonian matrix.
c            
            do 40 i=1,nend
               do 50 j=nbeg,nend
                  b(i,j) = cdotc(n,pvec(1,i),1,hpvec(1,j),1)
                  b(j,i) = cdotc(n,pvec(1,j),1,hpvec(1,i),1)
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
      return
 1    format(/,1x,'davidson eigenvalue solver using preconditioning',
     1                                                  /,10x,
     2            'number of roots               = ',i4,/,10x,
     3            'number of trials              = ',i4,/,10x,
     4            'maximum number of iterations  = ',i4,/,10x,
     5            'maximum number of vectors     = ',i4,/,10x,     
     6            'convergence criterion for rms = ',e15.8)      
 2    format(/,1x,'davidson eigenvalue solver',/,10x,
     1            'number of roots               = ',i4,/,10x,
     2            'number of trials              = ',i4,/,10x,
     3            'maximum number of iterations  = ',i4,/,10x,
     4            'maximum number of vectors     = ',i4,/,10x,
     5            'convergence criterion for rms = ',e15.8)      
 3    format(/,5x,'hamiltonian scale factor = ',e15.8)
 4    format(/,5x,'cannot even begin davidson calculation:',/,5x,
     1            'orthonormalization of initial vectors yields null'
     2            '  set')
 5    format(/,5x,'beginning davidson iterations:',/,5x,
     1            'error initialized at ',e15.8)
 6    format(/,1x,'beginning next davidson cycle',/,5x,
     1             'iteration = ',i4,1x,'size of vector space = ',i5) 
 7    format(/,5x,a80)
 8    format(/,30x,'summary',/,17x,'root',31x,'energy',/,46x,
     1                             'real',12x,'imag' )
 9    format(15x,i4,20x,e15.8,1x,e15.8)
 11   format(/,5x,'number of added vectors will exceed maxvec',/,5x,
     1            'number of vectors actually added = ',i4)           
 12   format(/1x,'no more orthonormal vectors can be added')
      end       
