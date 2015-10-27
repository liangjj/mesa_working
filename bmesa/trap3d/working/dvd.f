*deck dvd.f
c***begin prologue     dvd
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
c***end prologue       dvd
      subroutine dvd(hbuf,ibuf,diag,trials,etrial,eig01,eig02,eig03,
     1               u01,u02,u03,pvec,hpvec,vec,b,btmp,resid,eig,
     2               t1,t2,scale,cnverg,thresh,n1,n2,n3,dim,n,nroots,
     3               ntrials,maxit,lenbuf,drctv,incore,ntot,prnt)
      implicit integer (a-z)
      real*8 hbuf, diag, trials, etrial, pvec, hpvec, vec
      real*8 b, btmp, eig, eig01, eig02, eig03, u01, u02, u03
      real*8 resid, t1, t2, error, cnverg, thresh, sdot, scale, tmp
      real*8 zero, one, nrzero
      logical prnt, logkey, all, drctv, incore
      character*16 status
      character*3 itoc
      character*1 it
      character*80 title
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n), trials(n,ntrials)
      dimension etrial(ntrials), pvec(n,maxit), hpvec(n,maxit)
      dimension vec(n,maxit), b(maxit,maxit), btmp(maxit,maxit)
      dimension eig(ntrials), eig01(n1), eig02(n2), eig03(n3)
      dimension u01(n1,n1), u02(n2,n2), u03(n3,n3)
      dimension resid(n,maxit), t1(*), t2(*)
      dimension prnt(11)
      common/io/inp, iout
      data zero, one/ 0.d0, 1.d0/
      data nrzero / 1.0d-06 /

      it=itoc(dim)
      if(.not.incore) then
         call iosys('read integer "'//it//'d number of elements" '//
     1              'from ham',1,ntot,0,' ')
      endif
      if(drctv) then
         write(iout,1) nroots, ntrials, maxit, cnverg
      else
         write(iout,2) nroots, ntrials, maxit, cnverg
      endif
      if(prnt(1)) then
         title='trials'
         call prntrm(title,trials,n,ntrials,n,ntrials,iout)
      endif
      call copy(trials,pvec,n*ntrials)
      write(iout,3) scale
      nbeg=1
      nend=ntrials
      do 10 iter=1,maxit
c      
c        orthonormalize the new trials to the old vectors
c       
         call gschmt(pvec,thresh,n,nbeg,nend,nout,.true.)
         if(nout.ne.0) then
            nend=nbeg+nout-1
            if(prnt(11)) then
               write(iout,4) iter, nend
            endif
            if(prnt(3)) then
               call tstovl(pvec,n,nend)
            endif
            if(prnt(4)) then
               title='vectors iteration = '//itoc(iter)
               call prntrm(title,pvec(1,nbeg),n,nout,n,nout,iout)
            endif
c        operate with hamiltonian on these vectors
c
            call honv(hbuf,ibuf,diag,pvec(1,nbeg),hpvec(1,nbeg),n,
     1                nout,lenbuf,ntot,it,incore)
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
            if(prnt(11)) then
               do 60 i=1,nroots
                  tmp=scale*eig(i)
                  write(iout,6) i, tmp
 60            continue   
            endif
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
               if(prnt(11)) then
                  write(iout,5) i, error, status
               endif
 500        continue
            if(ncon.eq.nroots) then
               write(iout,8) nend
               do 2000 i=1,nroots
                  error = scale*sqrt (sdot(n,resid(1,i),1,
     1                                       resid(1,i),1) )
                  status='unconverged'
                  if(error.le.cnverg) then
                     status='converged'
                  endif
                  write(iout,9) i, scale*eig(i), error, status
 2000          continue   
               call iosys('write integer "size of davidson vector '//
     1                    'space" to ham',1,nend,0,' ')
               call vscale(eig,eig,scale,nend)
               call iosys('write real "davidson eigenvalues" to ham',
     1                     nend,eig,0,' ')
               call vscale(eig,eig,1.d0/scale,nend)
               call iosys('write real "davidson vectors" to ham',
     1                     nend*n,pvec,0,' ')
               return
            endif
            if(addvec.eq.0) then
               write(iout,7)
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
               call newvec(eig,diag,resid,pvec(1,nend+1),eig01,eig02,
     1                     eig03,u01,u02,u03,t1,t2,dim,n1,n2,n3,n,
     2                     addvec,drctv)
               nend=nend+addvec
               if(prnt(10)) then
                  title='new trial vectors iteration = '//itoc(iter)
                  call prntrm(title,pvec(1,nbeg),n,addvec,n,maxit,iout)
               endif
            endif
      else
         write(iout,11)
         do 3000 i=1,nroots
            error = scale*sqrt (sdot(n,resid(1,i),1,
     1                                 resid(1,i),1) )
            status='unconverged'
            if(error.le.cnverg) then
               status='converged'
            endif
            write(iout,9) i, scale*eig(i), error, status
 3000    continue   
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
     5            'convergence criterion for rms = ',e15.8)      
 2    format(/,1x,'davidson eigenvalue solver',/,10x,
     1            'number of roots               = ',i4,/,10x,
     2            'number of trials              = ',i4,/,10x,
     3            'maximum number of iterations  = ',i4,/,10x,
     4            'convergence criterion for rms = ',e15.8)      
 3    format(/,5x,'hamiltonian scale factor = ',e15.8)
 4    format(/,1x,'cycle = ',i3,5x,'size of vector space = ',i3) 
 5    format(/,5x,'root = ',i3,2x,'rms error = ',e15.8,2x,
     1            'status = ',a16)
 6    format(/,5x,'root = ',i3,2x,'energy = ',e15.8)
 7    format(/1x,'cannot add any more vectors.  will quit')
 8    format(/,1x,'all roots are converged.',/,1x,
     1            'number of davidson vectors = ',i5)
 9    format(/,5x,'root = ',i3,2x,'energy = ',e15.8,2x,
     1            'rms error = ',e15.8,/,5x,'status = ',a16)
 11   format(/1x,'no more orthonormal vectors can be added',/,1x,
     1           'summary of unconverged results' )
      end       



