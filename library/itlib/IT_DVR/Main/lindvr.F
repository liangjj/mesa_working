*deck lindvr.f
c***begin prologue     lindvr
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           iterative eigenvalue
c***author             schneider, barry (nsf)
c***source             
c***purpose            iterative linear system solver specialized 
c***                   to time-DVR using the Davidson algorithm.
c***                   
c***description        solve the linear set of equations
c***                        A |Psi> = |B>
c***                   where A and |B> are matrices.
c***
c***                   m = nz*ny*nx*nt*nc*2
c***                   n = nz*ny*nx*nt
c***references         
c
c***routines called    
c***end prologue       lindvr
      subroutine lindvr(hx,hy,hz,ht,eigx,eigy,eigz,ux,uy,uz,
     1                  v,pvec,hpvec,h,htmp,b,btmp,rhs,trials,resid,
     2                  soln,t1,t2,lufac,luind,scale,cnverg,thresh,eps,
     3                  nx,ny,nz,nt,nc,dim,m,n,nrhs,ntrial,maxit,
     4                  maxvec,prnt,precon)
      implicit integer (a-z)
      real*8 hx, hy, hz, ht, eigx, eigy, eigz, ux, uy, uz, v
      real*8 pvec, hpvec, h, htmp, b, btmp, rhs, trials
      real*8 resid, soln, t1, t2, error, cnverg, thresh, sdot, scale
      real*8 maxerr, eps
      real*4 time(20), delta(20), secnds
      complex*16 lufac
      logical prnt
      character*5 itoc
      character*1 it
      character*80 title
      character*(*) precon
      dimension pvec(m,maxvec), hpvec(m,maxvec)
      dimension rhs(m,nrhs), trials(m,ntrial), v(*)
      dimension h(maxvec,maxvec), htmp(maxvec,maxvec)
      dimension b(maxvec,nrhs), btmp(maxvec,nrhs)
      dimension hx(nx,nx), hy(ny,ny), hz(nz,nz), ht(nt,nt)
      dimension eigx(nx), eigy(ny), eigz(nz)
      dimension ux(nx,nx), uy(ny,ny), uz(nz,nz)
      dimension resid(m,maxvec), t1(*), t2(*), lufac(*), luind(*)
      dimension soln(m,nrhs)
      dimension prnt(11)
      common/io/inp, iout
      write(iout,1) m, nrhs, maxit, maxvec, ntrial, cnverg, thresh, eps
      write(iout,2) scale
c-----------------------------------------------------------------------c
c                                                                       c
c                    Initialization Section                             c
c                                                                       c
c-----------------------------------------------------------------------c
c
c     initialize the set of input vectors as the 
c
      if(prnt(2)) then
         title='input trial vectors'
         call prntrm(title,trials,m,ntrial,m,ntrial,iout)
         title='right hand side'
         call prntrm(title,rhs,m,nrhs,m,nrhs,iout)
         title='initial solution'
         call prntrm(title,soln,m,nrhs,m,nrhs,iout)
      endif
      call copy(trials,pvec,m*ntrial)
      time(1)=secnds(0.0)
      nbeg=1
      nend=ntrial
      call gschmt(pvec,thresh,m,1,nend,nout,.true.,.false.)
      time(2)=secnds(0.0)
      delta(1)=time(2)-time(1)      
      if(nend.eq.0) then
         call dvderr(1)
      else
         nend=nout
         if(prnt(2)) then
            title='initial trial vectors'
            call prntfm(title,pvec(1,nbeg),m,nout,m,nout,iout)
         endif
      endif
c    initialize the effect of the hamiltonian on these vectors.
c
      title='h on initial vectors'
      call htonv(hx,hy,hz,ht,v,pvec(1,nbeg),hpvec(1,nbeg),n,nc,
     1           nx,ny,nz,nt,nout,dim)
      time(3)=secnds(0.0)
      delta(2)=time(3)-time(2)
c
c     initialize the small hamiltonian matrix and right hand side.
c            
      call hinit(h,htmp,b,btmp,pvec,hpvec,rhs,m,nrhs,nend,maxvec)
      time(4)=secnds(0.0)
      delta(3)=time(4)-time(3)
      if(prnt(4)) then
         title='initial small matrix'
         call prntfm(title,htmp,nend,nend,maxvec,maxvec,iout)
         title='initial small right hand side'
         call prntfm(title,btmp,nend,nrhs,maxvec,maxvec,iout)
      endif
      write(iout,7) (delta(i),i=1,3)
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
      do while ( error.gt.cnverg.and.iter.lt.maxit )
         time(5)=secnds(0.0)
         iter = iter + 1
         write(iout,4) iter, nend
c     Step 1:
c
c           get solutions of the small matrix.
c
c                btmp holds the initial matrix which is destroyed.
c                note that resid is used as temporary storage in vscale.
c
         title='iteration = '//itoc(iter)//' solving linear system '
     1                      //'of dimension N = '//itoc(nend)
         write(iout,5) title
c        
         call lslv(htmp,btmp,resid,nend,nrhs,maxvec)
         call newsol(soln,pvec,btmp,m,nend,nrhs,maxvec)
         if(prnt(5)) then
            title='solutions of small matrix iteration = '//itoc(iter) 
            call prntfm(title,btmp,nend,nrhs,maxvec,nrhs,iout)
         endif
         time(6)=secnds(0.0)
         delta(4)=time(6)-time(5) 
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
         call lres(pvec,hpvec,btmp,rhs,scale,cnverg,resid,
     1             maxerr,soln,t1,m,nend,nrhs,con,uncon,
     2             maxvec,iter,prnt(5))
         time(7)=secnds(0.0)
         delta(5)=time(7)-time(6)
         if(maxerr.le.eps) then
            write(iout,*) 'writing guess vectors to ham'
            call iosys('write real "guess vectors" to ham',
     1                  nrhs*m,t1,0,' ') 
         endif
         if(con.eq.nrhs) then
c            call copy(t1,soln,m*nrhs)            
            call copy(soln,rhs,m*nrhs)
c
c           all solutions are converged.  copy them in to rhs and quit
c
            title='final solution'
            call prntrm(title,rhs,m,nrhs,m,nrhs,iout)
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
            write(iout,6) addvec
            chkno = nend + addvec
            if(chkno.le.maxvec) then
c
c              maximum number of vectors is still within the allowed
c              limits.  add vectors to the set from the unconverged
c              residuals and put them after the current vectors.
c
               nbeg = nend + 1
               if(precon.eq.'none') then
                  if(dim.eq.1) then
                     call v1t2h0(resid,t1,ux,nx,nt,nc,addvec)
                     call v1te(t1,eigx,nx,nt,nc,addvec)
                     call v1t2dvr(t1,pvec(1,nbeg),ux,nx,nt,nc,addvec)
                  elseif(dim.eq.2) then
                     call v2t2h0(resid,t1,ux,uy,nx,ny,nt,nc,addvec)
                     call v2te(t1,eigx,eigy,nx,ny,nt,nc,addvec)
                     call v2t2dvr(t1,pvec(1,nbeg),ux,uy,nx,ny,nt,
     1                            nc,addvec)
                  elseif(dim.eq.3) then
                     call v3t2h0(resid,t1,ux,uy,uz,nx,ny,nz,nt,
     1                           nc,addvec)
                     call v3te(t1,eigx,eigy,eigz,nx,ny,nz,nt,nc,addvec)
                     call v3t2dvr(t1,pvec(1,nbeg),ux,uy,uz,nx,ny,nz,
     1                                           nt,nc,addvec)
                  endif
               else
                  call cblslv(resid,pvec(1,nbeg),t1,lufac,luind,n,nc,
     1                        addvec)
               endif
               time(8)=secnds(0.0)
               delta(6)=time(8)-time(7)
               nend=nend+addvec
c
c              orthonormalize the new trials to the old vectors
c              to get an additional nout vectors.
c        
               call gschmt(pvec,thresh,m,nbeg,nend,nout,.true.,.false.)
               time(9)=secnds(0.0)
               delta(7)=time(9)-time(8)
c               
               if(nout.eq.0) then
c
c                 no more vectors write out unconverged results.
c               
c                  call dvderr(2)
                  return
               endif   
               nend=nbeg+nout-1
c               if(prnt(10)) then
c                  call tstovl(pvec,m,nend)
c               endif
c
c           we now have either a reconstituted set of vectors or
c           an augmented set of vectors.               
c
c              operate with hamiltonian on these vectors
c
               title='h on initial vectors'
               call htonv(hx,hy,hz,ht,v,pvec(1,nbeg),hpvec(1,nbeg),
     1                   n,nc,nx,ny,nz,nt,nout,dim)
               time(10)=secnds(0.0)
               delta(8)=time(10)-time(9)
c
c           update the curent small hamiltonian matrix and right hand sides.
c            
               call hupdat(h,htmp,b,btmp,pvec,hpvec,rhs,m,nrhs,nbeg,
     1                     nend,maxvec)
               time(11)=secnds(0.0)
               delta(9)=time(11)-time(10)
            else
               call dvderr(3)
            endif
         endif
         write(iout,8) iter, (delta(i),i=4,9)            
      enddo        
      call iosys('write integer "size of davidson vector '//
     1           'space" to ham',1,nend,0,' ')
      call iosys('write real "davidson vectors" to ham',
     1            nend*m,pvec,0,' ')            
      return
 1    format(/,1x,'davidson linear system solver using preconditioning',
     1                                                  /,10x,
     2            'number of equations                 = ',i4,/,10x, 
     3            'number of right hand sides          = ',i4,/,10x,
     4            'maximum number of iterations        = ',i4,/,10x,
     5            'maximum number of vectors           = ',i4,/,10x,     
     6            'number of initial trial vectors     = ',i4,/,10x,     
     7            'convergence criterion for rms       = ',e15.8,/,10x,
     8            'overlap threshold criterion         = ',e15.8,/,10x,      
     9            'restart criterion                   = ',e15.8)      
 2    format(/,5x,'hamiltonian scale factor = ',e15.8)
 3    format(/,5x,'beginning davidson iterations:',/,5x,
     1            'error initialized at ',e15.8) 
 4    format(/,1x,'beginning next davidson cycle',/,5x,
     1             'iteration = ',i4,1x,'size of vector space = ',i5)
 5    format(/,5x,a80)
 6    format(/,5x,'number of vectors actually added = ',i4)
 7    format(/,1x,'time for schmidt orthogonalization of initial '
     1            'vectors = ',e15.8,
     2       /,1x,'time for H on initial vectors          '
     3            '               = ',e15.8,
     4       /,1x,'time for initialization of small  hamiltonian       '
     5            '  = ',e15.8)
 8    format(/,1x,'timing at iteration = ',i4,
     1       /,1x,'time for solution of small set of linear '
     2            'equations = ',e15.8,
     3       /,1x,'time for residual calculation                      '
     4            '= ',e15.8,
     5       /,1x,'time for preconditioning                         '
     6            '  = ',e15.8,
     7       /,1x,'time for schmidt orthogonalization              '
     8            '   = ',e15.8,
     9       /,1x,'time for H on vectors                           '
     x            '   = ',e15.8,
     x       /,1x,'time to update small matrix                     '
     x            '   = ',e15.8)
      end       

