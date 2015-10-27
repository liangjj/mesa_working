*deck dvd1d.f
c***begin prologue     dvd1d
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
c***end prologue       dvd1d
      subroutine dvd1d(ham,diag,v,trials,pvec,hpvec,vec,b,
     1                  btmp,eig,eig0,u,resid,scale,cnverg,thresh,
     2                  n,nroots,ntrials,maxit,n0,drctv,ops)
      implicit integer (a-z)
      real*8 ham, diag, trials, v, pvec, hpvec, vec
      real*8 b, btmp, eig, eig0, u, resid, error, cnverg
      real*8 thresh, sdot, test, scale, tmp
      real*8 zero, one, nrzero
      logical prnt, logkey, all, drctv
      character*(*) ops
      character*16 status
      character*3 itoc
      character*80 title
      dimension ham(n,n), diag(n), v(n), trials(n,ntrials)
      dimension pvec(n,maxit), hpvec(n,maxit)
      dimension vec(n,maxit), b(maxit,maxit), btmp(maxit,maxit)
      dimension eig(maxit), eig0(n0), u(n,n0), resid(n,maxit)
      dimension prnt(10)
      common/io/inp, iout
      data zero, one/ 0.d0, 1.d0/
      data nrzero / 1.0d-06 /
      prnt(1)=logkey(ops,'print=eigen=trials',.false.,' ')
      prnt(2)=logkey(ops,'print=eigen=hamiltonian',.false.,' ')
      prnt(3)=logkey(ops,'print=eigen=overlaps',.false.,' ')
      prnt(4)=logkey(ops,'print=eigen=vectors',.false.,' ')
      prnt(5)=logkey(ops,'print=eigen=h-on-vectors',.false.,' ')
      prnt(6)=logkey(ops,'print=eigen=small-matrix',.false.,' ')
      prnt(7)=logkey(ops,'print=eigen=transformed-vectors',.false.,' ')
      prnt(8)=logkey(ops,'print=eigen=transformed-h-on-vectors',
     1               .false.,' ')
      prnt(9)=logkey(ops,'print=eigen=residuals',.false.,' ')
      prnt(10)=logkey(ops,'print=eigen=new-trial-vectors',.false.,' ')
      all=logkey(ops,'print=eigen=all',.false.,' ')
      if(drctv) then
         write(iout,1) nroots, ntrials, maxit, cnverg, n0
      else
         write(iout,2) nroots, ntrials, maxit, cnverg
      endif
      if(all) then
         do 1500 i=1,10
            prnt(i)=.true.
 1500    continue   
      endif       
      if(prnt(1)) then
         title='trials'
         call prntrm(title,trials,n,ntrials,n,ntrials,iout)
      endif
      call copy(trials,pvec,n*ntrials)
      call preph(ham,diag,v,n)
      if(prnt(2)) then
         title='hamiltonian'
         call prntrm(title,ham,n,n,n,n,iout)
         title='diagonals'
         call prntrm(title,diag,n,1,,n,1,iout)
      endif
      write(iout,*) 'scale =',scale
      nbeg=1
      nend=ntrials
      do 10 iter=1,maxit
c      
c        orthonormalize the new trials to the old vectors
c       
         call gschmt(pvec,thresh,n,nbeg,nend,nout,.true.)
         if(nout.ne.0) then
            nend=nbeg+nout-1
            write(iout,3) iter, nend
            if(prnt(3)) then
               call tstovl(pvec,n,nend)
            endif
            if(prnt(4)) then
               title='vectors iteration = '//itoc(iter)
               call prntrm(title,pvec(1,nbeg),n,nout,n,nout,iout)
            endif
c        operate with hamiltonian on these vectors
c
            call honv(ham,diag,pvec(1,nbeg),hpvec(1,nbeg),n,nout)
            if(prnt(5)) then
               title='h on vectors iteration = '//itoc(iter)
               call prntrm(title,hpvec(1,nbeg),n,nout,n,nout,iout)
            endif
c
c        update the curent small hamiltonian matrix.
c            
            do 20 i=1,nend
               do 30 j=nbeg,nend
                  b(i,j) = sdot(n,pvec(1,i),1,hpvec(1,j),1)
                  b(j,i) = b(i,j)
 30            continue   
 20         continue
            do 40 i=1,nend   
               do 50 j=1,i
                  btmp(i,j) = b(i,j)
                  btmp(j,i) = b(i,j)
 50            continue
 40         continue  
            if(prnt(6)) then
               title='small matrix iteration = '//itoc(iter)
               call prntrm(title,btmp,nend,nend,maxit,maxit,iout)
            endif
c
c        diagonalize small matrix
c
            call tred2(maxit,nend,btmp,eig,resid,btmp)
            call tql2(maxit,nend,eig,resid,btmp,ierr)
            if(ierr.ne.0) then
               call lnkerr('error from direct diagonalization routine')
            endif
            do 60 i=1,nroots
               tmp=scale*eig(i)
               write(iout,5) i, tmp
 60         continue   
c                
c
c        transform vectors and hamiltonian on vectors to new basis
c
            call ebcxx(vec,pvec,btmp,n,nend,nend,n,n,maxit)
            if(prnt(7)) then
               title='transformed vectors iteration = '//itoc(iter)
               call prntrm(title,vec,n,nend,n,maxit,iout)
            endif
c
c        replace the old vectors by the new ones
c
            call copy(vec,pvec,n*nend)                 
c
c        calculate new hamiltonian on vectors
c
            call ebcxx(vec,hpvec,btmp,n,nend,nend,n,n,maxit)
            if(prnt(8)) then
               title='transformed  h on vectors iteration = '
     1                //itoc(iter)
               call prntrm(title,vec,n,nend,n,maxit,iout)
            endif
c
c        replace old hamiltonian on vectors by the new ones
c
            call copy(vec,hpvec,n*nend)                 
c   
c        calculate residuals for desired roots and test for convergence
c
            addvec=0
            ncon=0
            do 500 i=1,nroots
               do 600 j=1,n
                  resid(j,i) = hpvec(j,i) - eig(i)*pvec(j,i) 
 600           continue
               error = scale*sqrt (sdot(n,resid(1,i),1,resid(1,i),1) )
               if(error.le.cnverg) then
                  status='converged'
                  ncon=ncon+1
               else
                  status='unconverged'
                  addvec=addvec+1
                  call copy(resid(1,i),resid(1,addvec),n)
                  eig(addvec) = eig(i)
               endif
               write(iout,4) i, error, status
 500        continue
            if(ncon.eq.nroots) then
               write(iout,7)
               return
            endif
            if(addvec.eq.0) then
               write(iout,6)
               return
            else               
               if(prnt(9)) then
                  title='residuals iteration = '//itoc(iter)
                  call prntrm(title,resid,n,addvec,n,maxit,iout)
               endif
c
c        make sure we will not exceed maxit
c
               newnum = nend + addvec
               newnum = min(newnum,maxit)
               addvec = newnum - nend
c
c        reform hamiltonian in new basis
c
               do 700 i=1,nend
                  do 800 j=1,i
                     b(i,j) = sdot(n,pvec(1,i),1,hpvec(1,j),1)
                     b(j,i) = b(i,j)
 800              continue
 700           continue
               if(prnt(6)) then
                  title='reformed small matrix iteration = '//itoc(iter)
                  call prntrm(title,b,nend,nend,maxit,maxit,iout)
               endif
c
c            prepare new trial vectors from residuals for unconverged
c            roots
c
               nbeg = nend + 1         
               call newvec(eig,diag,resid,pvec(1,nend+1),eig0,
     1                     u,n,addvec,n0,drctv)
               nend=nend+addvec
               if(prnt(10)) then
                  title='new trial vectors iteration = '//itoc(iter)
                  call prntrm(title,pvec(1,nbeg),n,addvec,n,maxit,iout)
               endif
            endif
      else
         write(iout,6) 
         return
      endif
c           
 10   continue
      return
 1    format(/,1x,'davidson eigenvalue solver using preconditioning',
     1                                                  /,10x,
     2            'number of roots               = ',i4,/,10x,
     3            'number of trials              = ',i4,/,10x,
     4            'maximum number of iterations  = ',i4,/,10x,
     5            'convergence criterion for rms = ',e15.8,/,10x,
     6            'number of zero order vectors  = ',i4)      
 2    format(/,1x,'davidson eigenvalue solver',/,10x,
     1            'number of roots               = ',i4,/,10x,
     2            'number of trials              = ',i4,/,10x,
     3            'maximum number of iterations  = ',i4,/,10x,
     4            'convergence criterion for rms = ',e15.8)      

 3    format(/,1x,'cycle = ',i3,5x,'size of vector space = ',i3) 
 4    format(/,5x,'root = ',i3,2x,'rms error = ',e15.8,2x,
     1            'status = ',a16)
 5    format(/,5x,'root = ',i3,2x,'energy = ',e15.8)
 6    format(/1x,'no more vectors.  will quit')
 7    format(/,5x,'all roots are converged. quit')
 8    format(/,5x,'root = ',i3,2x,'energy = ',e15.8,2x,'rms error = ',
     1                                        e15.8,2x,'status = ',a16)
      end       



