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
      subroutine nldvd(hbuf,ibuf,diag,vnl,trials,etrial,eig01,eig02,
     1                 eig03,u01,u02,u03,pvec,hpvec,vec,b,btmp,resid,
     2                 eig,t1,t2,scale,cnverg,thresh,n1,n2,n3,dim,n,
     3                 nroots,ntrials,maxit,lenbuf,nel,drctv,
     4                 incore)
      implicit integer (a-z)
      real*8 hbuf, diag, vnl, trials, etrial, pvec, hpvec, vec
      real*8 resid, eig01, eig02, eig03, u01, u02, u03
      real*8 eig, b, btmp, t1, t2
      real*8 scale, cnverg, thresh, sdot
      real*8 error, tmp
      character*80 title
      character*3 itoc
      character*1 it
      character*16 status
      logical drctv, prnt
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n), trials(n,ntrials)
      dimension etrial(ntrials), pvec(n,maxit), hpvec(n,maxit)
      dimension vec(n,maxit), resid(n,maxit), vnl(n)
      dimension b(maxit,maxit), btmp(maxit,maxit)
      dimension eig(n)
      dimension p1(n1,0:n1-1), p2(n2,0:n2-1)
      dimension p3(n3,0:n3-1)
      dimension eig01(n1), eig02(n2), eig03(n3)
      dimension u01(n1,n1), u02(n2,n2), u03(n3,n3)
      dimension t1(*), t2(*)
      dimension prnt(11)
      common/io/inp, iout
      data pi/3.14159265358979323844d0/
      data zero, one/ 0.d0, 1.d0/
      data nrzero / 1.0d-06 /
      data first / 0 /
      save first, prnt, it
      if(first.eq.0) then
         it=itoc(dim)
         call iosys('read integer "davidson print options" from ham',
     1               11,prnt,0,' ')
         if(drctv) then
            write(iout,1) nroots, ntrials, maxit, cnverg
         else
            write(iout,2) nroots, ntrials, maxit, cnverg
         endif
         first = 1
      endif
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
            if(prnt(11)) then
               write(iout,6) iter, nend
            endif
            if(prnt(3)) then
               call tstovl(pvec,n,nend)
            endif
            if(prnt(4)) then
               title='vectors iteration = '//itoc(iter)
               call prntfm(title,pvec(1,nbeg),n,nout,n,nout,iout)
            endif
            call honv(hbuf,ibuf,diag,pvec(1,nbeg),
     1                hpvec(1,nbeg),n,nout,lenbuf,nel,
     2                it,incore)
            call addvnl(pvec(1,nbeg),hpvec(1,nbeg),vnl,n,nout)
            if(prnt(5)) then
               title='h on vectors iteration = '//itoc(iter)
               call prntfm(title,hpvec(1,nbeg),n,nout,n,nout,iout)
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
               call prntfm(title,btmp,nend,nend,maxit,maxit,iout)
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
                  tmp=eig(i)/scale
                  write(iout,7) i, tmp
 60            continue   
            endif
c                
c
c        transform vectors and hamiltonian on vectors to new basis
c
            call ebcxx(vec,pvec,btmp,n,nend,nend,n,n,maxit)
            if(prnt(7)) then
               title='transformed vectors iteration = '//itoc(iter)
               call prntfm(title,vec,n,nend,n,maxit,iout)
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
               call prntfm(title,vec,n,nend,n,maxit,iout)
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
               error=error/scale
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
                  write(iout,8) i, error, status
               endif
 500        continue
            if(ncon.eq.nroots) then
               write(iout,9) nend
               do 2000 i=1,nroots
                  error = sqrt (sdot(n,resid(1,i),1,
     1                                 resid(1,i),1) )
                  status='unconverged'
                  if(error.le.cnverg) then
                     status='converged'
                  endif
                  write(iout,12) i, eig(i)/scale, error/scale, status
 2000          continue   
               return
            endif
            if(addvec.eq.0) then
               write(iout,11)
               return
            else               
               if(prnt(9)) then
                  title='residuals iteration = '//itoc(iter)
                  call prntfm(title,resid,n,addvec,n,maxit,iout)
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
                  call prntfm(title,b,nend,nend,maxit,maxit,iout)
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
                  call prntfm(title,pvec(1,nbeg),n,addvec,n,maxit,iout)
               endif
            endif
      else
         write(iout,11) 
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
 6    format(/,1x,'davidson cycle = ',i3,5x,'size of vector space = ',
     1                                i3) 
 7    format(/,5x,'root = ',i3,2x,'davidson energy = ',e15.8)
 8    format(/,5x,'root = ',i3,2x,'root rms error = ',e15.8,2x,
     1            'status = ',a16)
 9    format(/,1x,'all roots are converged.',/,1x,
     1            'number of davidson vectors = ',i5)
 12   format(/,10x,'root = ',i3,2x,'energy = ',e15.8,2x,
     1            'rms error = ',e15.8,/,10x,'status = ',a16)
 11   format(/1x,'no more vectors.  will quit')
      end       
