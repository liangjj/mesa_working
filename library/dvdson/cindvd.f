*deck cindvd.f
c***begin prologue     cindvd
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative eigenvalue
c***author             schneider, barry (nsf)
c***source             
c***purpose            iterative linear system solver using 
c***                   Davidson algorithm for complex matrices.
c***                   
c***description        solve the linear set of equations
c***                   ( E - H ) |Psi> = |B>
c***                   where |B> is a matrix.
c***references         
c
c***routines called    
c***end prologue       lindvd
      subroutine cindvd(hbuf,ibuf,diag,energy,eig01,eig02,eig03,
     1                  u0l1,u0l2,u0l3,u0r1,u0r2,u0r3,pvec,hpvec,vec,
     2                  rhs,h,htmp,b,btmp,resid,t1,t2,
     3                  scale,cnverg,thresh,n1,n2,n3,dim,n,nrhs,maxit,
     4                  maxvec,lenbuf,drctv,incore,ntot,prnt,mattyp)
      implicit integer (a-z)
      complex*16 hbuf, diag, pvec, hpvec, vec, rhs
      complex*16 h, htmp, b, btmp, eig01, eig02, eig03
      complex*16 u0r1, u0r2, u0r3, u0l1, u0l2, u0l3
      complex*16 resid, t1, t2, cdotc, tmp
      real*8 energy, error, cnverg, thresh, scale, maxerr
      logical prnt, drctv, incore
      character*(*) mattyp
      character*5 itoc
      character*1 it
      character*80 title
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n)
      dimension pvec(n,maxvec), hpvec(n,maxvec), rhs(n,nrhs)
      dimension vec(n,maxvec), h(maxvec,maxvec), htmp(maxvec,maxvec)
      dimension b(maxvec,nrhs), btmp(maxvec,nrhs)
      dimension eig01(n1), eig02(n2), eig03(N3)
      dimension u0l1(n1,n1), u0l2(n2,n2), u0l3(n3,n3)
      dimension u0r1(n1,n1), u0r2(n2,n2), u0r3(n3,n3)
      dimension resid(n,maxvec), t1(n,*), t2(n,*)
      dimension prnt(11)
      common/io/inp, iout
      it=itoc(dim)
      if(.not.incore) then
         call iosys('read integer "'//it//
     1              'd number of elements" from ham',1,ntot,0,' ')
      endif
      if(drctv) then
         write(iout,1) nrhs, maxit, maxvec, cnverg
      else
         write(iout,2) nrhs, maxit, maxvec, cnverg
      endif
      write(iout,3) scale
c-----------------------------------------------------------------------c
c                                                                       c
c                    Initialization Section                             c
c                                                                       c
c-----------------------------------------------------------------------c
c
c     initialize the set of input vectors as the 
c     orthonormalized set of right hand side vectors
c
      call cc2opy(rhs,pvec,n*nrhs)
      nbeg=1
      nend=nrhs
      call cschmt(pvec,thresh,n,1,nend,nout,.true.)
      if(nend.eq.0) then
         write(iout,4)
         call lnkerr('quit davidson')
      else
         nend=nout
         if(prnt(2)) then
            title='initial trial vectors'
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
c     initialize the small hamiltonian matrix and right hand side.
c            
      do 10 i=1,nend
         do 20 j=1,i
            h(i,j) = - cdotc(n,pvec(1,i),1,hpvec(1,j),1)
            h(j,i) = - cdotc(n,pvec(1,j),1,hpvec(1,i),1)
            htmp(i,j) = h(i,j)
            htmp(j,i) = h(i,j)
 20      continue   
         h(i,i) = energy + h(i,i)
         htmp(i,i) = h(i,i)
         do 30 j=1,nrhs
            b(i,j) = cdotc(n,pvec(1,i),1,rhs(1,j),1)
            btmp(i,j) = b(i,j)
 30      continue   
 10   continue
      if(prnt(4)) then
         title='initial small matrix'
         call prntcm(title,htmp,nend,nend,maxvec,maxvec,iout)
         title='initial small right hand side'
         call prntcm(title,btmp,nend,nrhs,maxvec,maxvec,iout)
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
      iter=0
      error=1.d+10
      write(iout,5) error
      do while ( error.gt.cnverg.and.iter.lt.maxit )
         iter = iter + 1
         write(iout,6) iter, nend
c     Step 1:
c
c           get solutions of the small matrix.
c
c                btmp holds the initial matrix which is destroyed.
c                note that resid is used as temporary storage in vscale.
c
         title='iteration = '//itoc(iter)//' solving linear system '
     1                      //'of dimension N = '//itoc(nend)
         write(iout,7) title
c        
         if(prnt(4)) then
            title='small matrix'
            call prntcm(title,htmp,nend,nend,maxvec,maxvec,iout)
            title='small right hand side'
            call prntcm(title,btmp,nend,nrhs,maxvec,nrhs,iout)
         endif
         call lsolve(htmp,htmp,btmp,btmp,resid,nend,nrhs,
     1               maxvec,mattyp)
         if(prnt(5)) then
            title='solutions of small matrix iteration = '//itoc(iter) 
            call prntcm(title,btmp,nend,nrhs,maxvec,nrhs,iout)
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
         call clares(pvec,hpvec,btmp,rhs,energy,scale,cnverg,resid,
     1               maxerr,t1,t2,n,nend,nrhs,con,uncon,
     2               maxvec,iter,prnt(5))
         if(con.eq.nrhs) then
c
c           all solutions are converged.  copy them in to rhs and quit
c
            call cc2opy(t1,rhs,n*con)            
            write(iout,8)
            return
         else
c
c           all solutions are not converged.  set the error to the largest
c           current error and continue the iteration sequence or 
c           quit if maxvec is exceeded.
c
            error=min(error,maxerr)
c           how many new vectors could be added in principle
c
            numnew = maxvec - nend
c                                             
c           how many will we add
c
            addvec = min(numnew,uncon)
            write(iout,9) addvec
            chkno = nend + addvec
            if(chkno.le.maxvec) then
c
c              maximum number of vectors is still within the allowed
c              limits.  add vectors to the set from the unconverged
c              residuals and put them after the current vectors.
c
               nbeg = nend + 1
               call clavec(energy,diag,resid,pvec(1,nbeg),eig01,eig02,
     1                     eig03,u0l1,u0l2,u0l3,u0r1,u0r2,u0r3,
     2                     t1,t2,dim,n1,n2,n3,n,addvec,drctv,
     3                     iter,prnt(9))
               nend=nend+addvec
c
c              orthonormalize the new trials to the old vectors
c              to get an additional nout vectors.
c        
               call cschmt(pvec,thresh,n,nbeg,nend,nout,.true.)
c               
               if(nout.eq.0) then
c
c                 no more vectors write out unconverged results.
c               
                  write(iout,11)
                  return
               endif   
               nend=nbeg+nout-1
               if(prnt(10)) then
                  call tstovl(pvec,pvec,n,nend,mattyp)
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
c           update the curent small hamiltonian matrix and right hand sides.
c            
               do 40 i=1,nend
                  do 50 j=nbeg,nend
                      h(i,j) = - cdotc(n,pvec(1,i),1,hpvec(1,j),1)
                      h(j,i) = - cdotc(n,pvec(1,j),1,hpvec(1,i),1)
 50               continue
                  do 60 j=1,nrhs
                     b(i,j) = cdotc(n,pvec(1,i),1,rhs(1,j),1)   
 60               continue   
 40            continue
               do 70 i=nbeg,nend
                  h(i,i) = energy + h(i,i)
 70            continue   
               do 80 i=1,nend   
                  do 90 j=1,i
                     htmp(i,j) = h(i,j)
                     htmp(j,i) = h(i,j)
 90               continue
                  do 100 j=1,nrhs
                     btmp(i,j) = b(i,j)
 100              continue   
 80            continue  
            else
               write(iout,12)
            endif
         endif            
      enddo        
      return
 1    format(/,1x,'davidson linear system solver using preconditioning',
     1                                                  /,10x,
     2            'number of right hand sides    = ',i4,/,10x,
     3            'maximum number of iterations  = ',i4,/,10x,
     4            'maximum number of vectors     = ',i4,/,10x,     
     5            'convergence criterion for rms = ',e15.8)      
 2    format(/,1x,'davidson linear system solver',/,10x,
     1            'number of right hand sides    = ',i4,/,10x,
     2            'maximum number of iterations  = ',i4,/,10x,
     3            'maximum number of vectors     = ',i4,/,10x,
     4            'convergence criterion for rms = ',e15.8)      
 3    format(/,5x,'hamiltonian scale factor = ',e15.8)
 4    format(/,5x,'cannot even begin davidson calculation:',/,5x,
     1            'orthonormalization of initial vectors yields null'
     2            '  set')
 5    format(/,5x,'beginning davidson iterations:',/,5x,
     1            'error initialized at ',e15.8) 
 6    format(/,1x,'beginning next davidson cycle',/,5x,
     1             'iteration = ',i4,1x,'size of vector space = ',i5)
 7    format(/,5x,a80)
 8    format(/,5x,'all solutions are converged. quit')
 9    format(/,5x,'number of added vectors will exceed maxvec',/,5x,
     1            'number of vectors actually added = ',i4)
 11   format(/,5x,'cannot continue since maxvec is exceeded')
 12   format(/1x,'no more orthonormal vectors can be added')
      end       
