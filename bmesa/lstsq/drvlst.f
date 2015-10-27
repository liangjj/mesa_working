*deck drvlst.f
c***begin prologue     drvlst
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           least squares
c***author             schneider, barry (nsf)
c***source             
c***purpose            least squares solution of linear equations.
c***                   
c***references         
c
c***routines called    
c***end prologue       drvlst
      subroutine drvlst(ham,diag,rhs,trials,pvec,hpvec,vec,tvec,b,
     1                  btmp,sol,soltmp,energy,cnvrg,thresh,
     2                  ipvt,n,nrhs,ntrials,maxit,ops)
      implicit integer (a-z)
      real*8 ham, diag, rhs, trials, pvec, hpvec, vec, tvec
      real*8 b, btmp, sol, soltmp, error, cnvrg
      real*8 thresh, sdot, maxerr, energy
      real*8 zero, one
      character*80 title
      character*(*) ops
      logical prnt, schmdt, all, logkey
      dimension ham(n,n), diag(n), rhs(n,nrhs), trials(n,ntrials)
      dimension pvec(n,maxit), hpvec(n,maxit)
      dimension vec(n,nrhs), tvec(n,nrhs), b(maxit,maxit)
      dimension btmp(maxit,maxit), sol(maxit,nrhs), soltmp(maxit,nrhs)
      dimension ipvt(maxit), prnt(7)
      common/io/inp, iout
      data zero, one/0.d0,1.d0/
      prnt(1)=logkey(ops,'print=linslv=trials',.false.,' ')
      prnt(2)=logkey(ops,'print=linslv=right-hand-sides',.false.,' ')
      prnt(3)=logkey(ops,'print=linslv=vectors',.false.,' ')
      prnt(4)=logkey(ops,'print=linslv=h-on-vectors',.false.,' ')
      prnt(5)=logkey(ops,'print=linslv=small-matrix',.false.,' ')
      prnt(6)=logkey(ops,'print=linslv=solutions',.false.,' ')
      prnt(7)=logkey(ops,'print=linslv=residuals',.false.,' ')
      all=logkey(ops,'print=linslv=all',.false.,' ')
      schmdt=logkey(ops,'linslv=two-schmidt-orthogonalizations',
     1              .false.,' ')
      if (all) then
          do 10 i=1,7
             prnt(i)=.true.
 10       continue
      endif 
      write(iout,1) nrhs, ntrials, maxit, cnvrg
c     setup the matrix and right hand side in convenient form.
c     after the call the diagonals are zeroed and the remaining
c     elements divided by the diagonals.  the same thing is done to the
c     right hand side.
      call copy(trials,pvec,n*ntrials)
      if (prnt(1)) then
          title='trials'
          call prntrm(title,pvec,n,ntrials,n,maxit,iout)
      endif
      call precon(ham,diag,rhs,pvec,energy,n,nrhs,ntrials)
      if (prnt(2)) then
          title='right hand sides'
          call prntrm(title,rhs,n,nrhs,n,nrhs,iout)
      endif
      nbeg=1
      nend=ntrials
      do 20 iter=1,maxit
         call gschmt(pvec,thresh,n,nbeg,nend,nout,schmdt)
c         call tstovl(pvec,n,nend)
         if(nout.ne.0) then
            nend=nbeg+nout-1
            if(prnt(3)) then
               title='new vectors'
               call prntrm(title,pvec(1,nbeg),n,nout,n,nout,iout)
            endif
c        operate with hamiltonian on this vector
c
         call honv(ham,pvec(1,nbeg),hpvec(1,nbeg),n,nout)
         if(prnt(4)) then
            title='hamiltonian on vectors'
            call prntrm(title,hpvec(1,nbeg),n,nout,n,nout,iout)
         endif
c
c        update the curent coefficient matrix and right hand side.
c            
         do 30 i=1,nbeg-1
            do 40 j=nbeg,nend
               b(i,j) = sdot(n,pvec(1,i),1,hpvec(1,j),1)
               b(j,i) = sdot(n,pvec(1,j),1,hpvec(1,i),1)
 40         continue   
 30      continue
         do 50 i=1,nrhs
            do 60 j=nbeg,nend
               sol(j,i)=sdot(n,pvec(1,j),1,rhs(1,i),1)            
 60         continue   
 50      continue   
         do 70 i=nbeg,nend
            do 80 j=nbeg,nend
               b(i,j) = sdot(n,pvec(1,i),1,hpvec(1,j),1)
 80         continue
 70      continue   
         do 90 i=nbeg,nend
            b(i,i) = b(i,i) + one
 90      continue   
c
         do 100 i=1,nrhs
            do 110 j=1,nend
               soltmp(j,i)=sol(j,i)
 110        continue   
 100      continue
         do 200 i=1,nend   
            do 210 j=1,nend
               btmp(i,j) = b(i,j)
 210        continue
 200     continue
         if(prnt(5)) then
            title='small matrix'
            call prntrm(title,b,nend,nend,maxit,maxit,iout)
            title='right hand side'
            call prntrm(title,sol,nend,nrhs,maxit,nrhs,iout)
         endif   
c
c        solve linear equations
c
         call sgefa(btmp,maxit,nend,ipvt,info)
         if(info.ne.0) then
            call lnkerr('matrix singular. quit')
         endif
         do 300 i=1,nrhs
            call sgesl(btmp,maxit,nend,ipvt,soltmp(1,i),0)
 300     continue   
c                
c
c        calculate solutions
c
         call ebcxx(vec,pvec,soltmp,n,nend,nrhs,n,n,maxit)
c
c        calculate hamiltonian on solution
c
         call ebcxx(tvec,hpvec,soltmp,n,nend,nrhs,n,n,maxit)
         if(prnt(6)) then
            title='solution matrix'
            call prntrm(title,vec,n,nrhs,n,maxit,iout)
         endif   
c
c        calculate residuals
c
         maxerr=0.d0
         do 400 i=1,nrhs
            do 410 j=1,n
               vec(j,i) = vec(j,i) + tvec(j,i) - rhs(j,i)
 410        continue
            error=sqrt ( sdot(n,vec(1,i),1,vec(1,i),1) )
            maxerr=max(maxerr,error)   
 400     continue
          if(prnt(7)) then
            title='residual matrix'
            call prntrm(title,vec,n,nrhs,n,nrhs,iout)
         endif
c
c        check convergence
c
         write(iout,2) iter, nend, maxerr
         if(maxerr.le.cnvrg) then
            write(iout,3)
            return
         endif  
c
c
c        prepare for next iteration
c
         call copy(hpvec(1,nbeg),pvec(1,nend+1),n*nout)
         nbeg = nend + 1
         nend = nbeg + nout - 1
      else
         write(iout,4) 
      endif
c           
 20   continue
      return
 1    format(/,1x,'least squares solution of linear equations',/,10x,
     1            'number of right hand sides    = ',i4,/,10x,
     2            'number of trials              = ',i4,/,10x,
     3            'maximum number of iterations  = ',i4,/,10x,
     4            'convergence criterion for rms = ',e15.8)
 2    format(/,5x,'iteration = ',i3,2x,'number of vectors = ',i3,2x,
     1            'max. rms error = ',e15.8)
 3    format(/,1x,'solution has converged to tolerance')
 4    format(/1x,'no more vectors.  will quit')      
      end       



