*deck qcklin.f
c***begin prologue    qcklin
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative eigenvalue
c***author             schneider, barry (nsf)
c***source             
c***purpose            iterative linear system solver specialized 
c***                   to time-DVR using the Davidson algorithm and
c***                   multigrid.
c***                   
c***description        solve the linear set of equations
c***                        A |Psi> = |B>
c***                   where A and |B> are matrices.
c***
c***                   m = nz*ny*nx*nt*2
c***                   ntot = nz*ny*nx*nt
c***                   n3d = nz*ny*nx
c
c***                   soln = guess input and final output solution
c***                   rhs  = right hand side
c
c***references         
c
c***routines called    
c***end prologue      qcklin
      subroutine qcklin(psoln,pdrver,pbasis,parray,pvoft,phamil,plufac,
     1                  pm,psln,ham,eig,trn,pvt,pdrv,scale,cnverg,
     2                  thresh,eps,maxit,maxvec,dim,spdim,nmax,msze,
     3                  ntot,n3d,mgr,grd,precon,nblck,qtitle,
     4                  prnt,mprnt)
      implicit integer (a-z)
      real*8 basis, array, voft, driver, soln, hamil, lufac
      real*8 scale, cnverg, thresh, eps, error, maxerr
      real*8 lin, rhs
      logical prnt, mprnt
      character*5 itoc
      character*1 it
      character*80 title
      character*(*) qtitle, precon
      dimension nmax(4,mgr), msze(mgr), ntot(mgr), n3d(mgr)
      dimension pm(mgr,mgr,dim)    
      dimension ham(mgr,dim), trn(mgr,dim), eig(mgr,dim)
      dimension pvt(mgr), pdrv(mgr), psln(mgr)
      dimension maxit(mgr), maxvec(mgr)
      dimension maxvec(*)< maxit(*), prnt(*), mprint(*)
      common/io/inp, iout
      pointer(pbasis,basis(1))
      pointer(parray,array(1))
      pointer(pvoft,voft(1))     
      pointer(pdrver,driver(1))
      pointer(psoln,soln(1))
      pointer(phamil,hamil(1)), (phamil,ihamil(1))
      pointer(plufac,lufac(1)), (plufac,indlu(1))
      pointer(prhs,rhs(1))
      pointer(plin,lin(1)), (plin,ilin(1))
c
c     define pointers to arrays for this grid and some needed integers
c
      nx=nmax(2,grd)
      ny=nmax(3,grd)
      nz=nmax(4,grd)
      nt=nmax(1,grd)
      n=n3d(grd)
      m=msze(grd)
      nsze=ntot(grd) 
      hx=ham(grd,2)
      hy=ham(grd,3)
      hz=ham(grd,4)
      ht=ham(grd,1)
      eigx=eig(grd,2)
      eigy=eig(grd,3)
      eigz=eig(grd,4)
      trnx=ham(grd,2)
      trny=ham(grd,3)
      trnz=ham(grd,4)
      v=pvt(grd)
      mxvc=maxvec(grd)
      mxit=maxit(grd)
      if(precon.ne.'none') then
         lumat=1
         luind=wpadti(lumat+2*nblck*nblck)
      endif
      h0=1
      h0pvt=wpadti(h0+msze(1)*msze(1))
c
c     get memory for this routine
c
c     davidson vectors = pvec
c
      pvec=1
c
c     davidson trials = trials
c
      trials=1
c
c     action of hamiltonian on pvec = hpvec
c
      hpvec=pvec+m*mxvc
c
c     hamiltonian projected into davidson space = h
c
      h=hpvec+m*mxvc
c
c     temporary storage of h = htmp so as not to recompute at each iteration
c
      htmp=h+mxvc*mxvc
c
c     right hand side projected into davidson space = b
c
      b=htmp+mxvc*mxvc
c
c     temporary storage of b = btmp so as not to recompute at each iteration
c
      btmp=b+mxvc
c
c     residual vector
c
      resid=btmp+mxvc
c
c     temporary storage
c
      t1=resid+m
c
      need=wpadti(t1+m)
c      
      call memory(need,plin,ngot,'lindvr',0)
c
      call copy(soln,lin(trials),m)      
      write(iout,1) qtitle
c-----------------------------------------------------------------------c
c                                                                       c
c                    Initialization Section                             c
c                                                                       c
c-----------------------------------------------------------------------c
      if(prnt(2)) then
         title='input trial vectors'
         call prntrm(title,lin(trials),m,1,m,1,iout)
         title='right hand side'
         call prntrm(title,rhs,m,1,m,1,iout)
         title='initial solution'
         call prntrm(title,soln,m,1,m,1,iout)
      endif
      call copy(lin(trials),lin(pvec),m)
      nbeg=1
      nend=nbeg
      call gschmt(lin(pvec),thresh,m,1,nend,nout,.true.)
      if(nend.eq.0) then
         call dvderr(1)
      else
         nend=nout
      endif
      if(prnt(2)) then
         title='initial trial vectors'
         call prntfm(title,lin(pvec),m,nout,m,nout,iout)
      endif
c
      call linit(array(ht),array(hx),array(hy),array(hz),voft(v),
     1           rhs,lin(pvec),lin(hpvec),lin(h),lin(htmp),lin(b),
     2           lin(btmp),n,m,nt,nx,ny,nz,nbeg,nend,nout,spdim,mxvc)
      if(prnt(4)) then
         title='initial small matrix'
         call prntfm(title,lin(htmp),nend,nend,mxvc,mxvc,iout)
         title='initial small right hand side'
         call prntfm(title,lin(btmp),nend,1,nend,1,iout)
      endif
c----------------------------------------------------------------------c
c                                                                      c
c                    Iteration Sequence                                c
c                                                                      c
c     iteration is continued until the solution is converged           c 
c     or if convergence is not achieved some preset maximum number of  c
c     iterations are performed.                                        c
c                                                                      c
c----------------------------------------------------------------------c
      error=1.d+10
      write(iout,2) error
      do 10 iter=1,mxit
         write(iout,3) iter, nend
c     Step 1:
c
c           get solutions of the small matrix.
c
c                btmp holds the initial matrix which is destroyed.
c                note that resid is used as temporary storage in vscale.
c
         title='iteration = '//itoc(iter)//' solving linear system '
     1                      //'of dimension N = '//itoc(nend)
         write(iout,1) title
c        
         call lslv(lin(htmp),lin(btmp),lin(resid),nend,1,maxvec)
         call newsol(soln,lin(pvec),lin(btmp),m,nend,1,mxvc)
c
         if(prnt(5)) then
            title='solutions of small matrix iteration = '//itoc(iter) 
            call prntfm(title,lin(btmp),nend,1,nend,1,iout)
         endif
c      
c     Step 2:
c
c           form the residuals and check for convergence.
c
c
c        unconverged residuals are in resid and the corresponding eigenvalues
c        are in etmp.
c 
         call lres(lin(pvec),lin(hpvec),lin(btmp),rhs,
     1             scale,cnverg,resid,maxerr,soln,lin(t1),
     2             m,nend,1,con,uncon,mxvc,iter,prnt(5))
         if(con.eq.1) then
c
c           all solutions are converged.
c
            write(iout,4) grd
            title='final solution'
            call prntrm(title,soln,m,1,m,1,iout)
            call memory(-ngot,plin,idum,'lindvr',idum)
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
            numnew = mxvc - nend
c                                             
c           how many will we add
c
            addvec = min(numnew,uncon)
            write(iout,5) addvec
            chkno = nend + addvec
            if(chkno.le.mxvc) then
c
c              maximum number of vectors is still within the allowed
c              limits.  add vectors to the set from the unconverged
c              residuals and put them after the current vectors.
c
               nbeg = nend + 1
               lpvec = pvec + (nbeg-1)*m
               lhpvec = hpvec + (nbeg-1)*m

c
c              get the next correction vector using the multigrid
c              procedure.
c
               if(precon.eq.'none') then
                  if(spdim.eq.1) then
                     call v12h0(lin(resid),lin(t1),array(trnx),
     1                          nx,nt,addvec)
                     call v1e(lin(t1),eig(eigx),nx,nt,addvec)
                     call v12dvr(lin(t1),lin(lpvec),array(trnx),
     1                           m,nx,nt,addvec)
                  elseif(spdim.eq.2) then
                     call v22h0(lin(resid),lin(t1),
     1                          array(trnx),array(trny),nx,ny,nt,addvec)
                     call v2e(lin(t1),eig(eigx),eig(eigy),
     1                        nx,ny,nt,addvec)
                     call v22dvr(lin(t1)lin(lpvec),
     1                           array(trnx),array(trny),array(trnz),
     2                           nx,ny,nt,addvec)
                  elseif(spdim.eq.3) then
                     call v32h0(lin(resid),lin(t1),array(trnx),
     1                          array(trny),array(trnz),
     2                          nx,ny,nz,nt,addvec)
                     call v3e(lin(t1),eig(eigx),eig(eigy),eig(eigz),
     1                        nx,ny,nz,nt,addvec)
                     call v32dvr(lin(t1),lin(lpvec),array(trnx),
     1                           array(trny),array(trnz),
     2                           nx,ny,nz,nt,addvec)
                  endif
               else
                  call blkslv(lin(resid),lin(lpvec),lin(t1),
     1                        lufac(luh),indlu(luind),nsze,addvec,grd)
               endif
               nend=nend+1
c
c              orthonormalize the new trials to the old vectors
c              to get an additional nout vectors.
c        
               call gschmt(lin(pvec),thresh,m,nbeg,nend,nout,.true.)
c               
               if(nout.eq.0) then
c
c                 no more vectors.  return with best you have.
c               
                  call memory(-ngot,plin,idum,'lindvr',idum)
                  return
               endif   
               nend=nbeg+nout-1
               if(prnt(10)) then
                  call tstovl(lin(pvec),lin(pvec),m,nend,mattyp)
               endif
c
c           we now have either a reconstituted set of vectors or
c           an augmented set of vectors.               
c
c              operate with hamiltonian on these vectors
c
               title='h on initial vectors'
               call honv(array(ht),array(hx),array(hy),array(hz),
     1                   voft(v),lin(lpvec),lin(lhpvec),n,nt,
     2                   nx,ny,nz,nout,spdim)
c
c           update the curent small hamiltonian matrix and right hand sides.
c            
               call hupdat(lin(h),lin(htmp),lin(b),lin(btmp),
     1                     lin(pvec),lin(hpvec),rhs,
     2                     m,1,nbeg,nend,maxvec)
            endif
         endif            
 10   continue   
      call memory(-ngot,plin,idum,'lindvr',idum)
      return
 1    format(/,5x,a80)
 2    format(/,5x,'beginning davidson iterations:',/,5x,
     1            'error initialized at ',e15.8) 
 3    format(/,5x,'beginning next davidson cycle',/,5x,
     1             'iteration = ',i4,1x,'size of vector space = ',i5)
 4    format(/,5x,'next vector obtained successfully on grid = ',i3)
 5    format(/,5x,'number of vectors actually added = ',i4)
      end       

