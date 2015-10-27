*deck nldvd.f
c***begin prologue     nldvd
c***date written       951229   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           non-linear schroedinger equation, gross-pitaevski,
c***                   davidson
c***author             schneider, barry (nsf)
c***source             timeprp
c***purpose            calculate eigenvalues of gross-pitaevski equation
c***                   using the davidson procedure
c***                   
c***references         
c
c***routines called    
c***end prologue       nldvd
      subroutine nldvd(hbuf,ibuf,diag,vnl,trials,etrial,pvec,hpvec,vec,
     1                 resid,psi,p1,p2,p3,eigr1,eigr2,eigr3,ham01,
     2                 ham02,ham03,eig01,eig02,eig03,u01,u02,u03,
     3                 b,btmp,eig,t1,t2,t3,table,natmin,natmax,hbar,
     4                 amass,omega,scatl,cnverg,thresh,nstep,number,
     5                 nd,dim,n,nroots,ntrials,maxit,lenbuf,drctv,
     6                 incore,prnt)
      implicit integer (a-z)
      real*8 hbuf, diag, vnl, trials, etrial, pvec, hpvec, vec, resid
      real*8 psi, p1, p2, p3, eigr1, eigr2, eigr3
      real*8 ham01, ham02, ham03, eig01, eig02, eig03
      real*8 b, btmp, eig, t1, t2, t3, table
      real*8 natmin, natmax, natoms
      real*8 hbar, amass, omega, scatl, cnverg, thresh
      real*8 gamma, norm0, pi, zero, one, delat
      real*8 crit, scale, u0, sdot, error, tmp
      character*80 title
      character*3 itoc
      character*1 it
      character*16 status
      logical drctv, prnt, all
      dimension nd(3)
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n), trials(n,ntrials)
      dimension etrial(ntrials), pvec(n,maxit), hpvec(n,maxit)
      dimension vec(n,maxit), resid(n,maxit), vnl(n)
      dimension psi(n), b(maxit,maxit), btmp(maxit,maxit)
      dimension eig(n)
      dimension p1(nd(1),0:nd(1)-1), p2(nd(2),0:nd(2)-1)
      dimension p3(nd(3),0:nd(3)-1)
      dimension eigr1(nd(1)), eigr2(nd(2)), eigr3(nd(3))
      dimension ham01(nd(1),nd(1)), ham02(nd(2),nd(2))
      dimension ham03(nd(3),nd(3))
      dimension eig01(nd(1)), eig02(nd(2)), eig03(nd(3))
      dimension u01(nd(1),nd(1)), u02(nd(2),nd(2)), u03(nd(3),nd(3))
      dimension t1(*), t2(*), t3(*)
      dimension table(number,3)
      dimension prnt(10)
      common/io/inp, iout
      data pi/3.14159265358979323844d0/
      data zero, one/ 0.d0, 1.d0/
      data nrzero / 1.0d-06 /
      it=itoc(dim)
      call iosys('read integer "'//it//'d number of elements" from ham',
     1            1,ntot,0,' ')
      norm0=1.d0/(4.d0*pi)
      norm0=sqrt(norm0)
      scale=omega*hbar
      u0=4.d0*pi*hbar*hbar*scatl/amass
      write(iout,1)
      write(iout,2) amass, omega, scatl
      if(drctv) then
         write(iout,3) nroots, ntrials, maxit, cnverg
      else
         write(iout,4) nroots, ntrials, maxit, cnverg
      endif
      write(iout,5) scale
      natoms=natmin
      delat=(natmax - natmin)/nstep
      count=0
      call copy(trials(1,1),psi,n)
      call smul(psi,psi,norm0,n)
      gamma=natoms*u0
      call potnl(vnl,eigr1,p1,psi,gamma,n)
      call copy(trials,pvec,n*ntrials)
      nbeg=1
      nend=ntrials
      do 10 iter=1,maxit
c      
c        orthonormalize the new trials to the old vectors
c       
         call gschmt(pvec,thresh,n,nbeg,nend,nout,.true.)
         if(nout.ne.0) then
            nend=nbeg+nout-1
            write(iout,6) iter, nend
            if(prnt(3)) then
               call tstovl(pvec,n,nend)
            endif
            if(prnt(4)) then
               title='vectors iteration = '//itoc(iter)
               call prntrm(title,pvec(1,nbeg),n,nout,n,nout,iout)
            endif
c
c        calculate non-linear term and then
c        operate with hamiltonian on these vectors
c            gamma=natoms*u0
c            call copy(pvec(1,1),psi,n)
c            call smul(psi,psi,norm0,n)
c            call potnl(vnl,eigr1,p1,psi,gamma,n)
            call hnlonv(hbuf,ibuf,diag,vnl,pvec(1,nbeg),
     1                  hpvec(1,nbeg),n,nout,lenbuf,ntot,
     2                  it,incore)
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
               tmp=eig(i)/scale
               write(iout,7) i, tmp
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
               error = sqrt (sdot(n,resid(1,i),1,resid(1,i),1) )
               error=error/(scale*norm0*norm0)
               if(error.le.cnverg) then
                  status='converged'
                  ncon=ncon+1
               else
                  status='unconverged'
                  addvec=addvec+1
                  call copy(resid(1,i),resid(1,addvec),n)
                  eig(addvec) = eig(i)
               endif
               write(iout,8) i, error, status
 500        continue
            if(ncon.eq.nroots) then
               write(iout,9)
               return
            endif
            if(addvec.eq.0) then
               write(iout,11)
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
     1                     eig03,u01,u02,u03,t1,t2,t3,dim,nd,n,
     2                     addvec,drctv)
               nend=nend+addvec
               if(prnt(10)) then
                  title='new trial vectors iteration = '//itoc(iter)
                  call prntrm(title,pvec(1,nbeg),n,addvec,n,maxit,iout)
               endif
            endif
      else
         write(iout,11) 
         return
      endif
c           
 10   continue
         natoms = natoms + delat
 1000 continue   
      return
 1    format(/,5x,'diagonalize time-independent non-linear schroedinger'
     1            ' equation')
 2    format(/,5x,'atomic mass        = ',e15.8,/,5x,
     1            'trap frequency     = ',e15.8,/,5x,
     2            'scattering length  = ',e15.8)
 3    format(/,1x,'davidson eigenvalue solver using preconditioning',
     1                                                  /,10x,
     2            'number of roots               = ',i4,/,10x,
     3            'number of trials              = ',i4,/,10x,
     4            'maximum number of iterations  = ',i4,/,10x,
     5            'convergence criterion for rms = ',e15.8)      
 4    format(/,1x,'davidson eigenvalue solver',/,10x,
     1            'number of roots               = ',i4,/,10x,
     2            'number of trials              = ',i4,/,10x,
     3            'maximum number of iterations  = ',i4,/,10x,
     4            'convergence criterion for rms = ',e15.8)      
 5    format(/,5x,'hamiltonian scale factor = ',e15.8)
 6    format(/,1x,'cycle = ',i3,5x,'size of vector space = ',i3) 
 7    format(/,5x,'root = ',i3,2x,'energy = ',e15.8)
 8    format(/,5x,'root = ',i3,2x,'rms error = ',e15.8,2x,
     1            'status = ',a16)
 9    format(/,5x,'all roots are converged. quit')
 11   format(/1x,'no more vectors.  will quit')
      end       
