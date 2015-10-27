*deck trmat.f
c***begin prologue     trmat
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           preconditioner
c***author             schneider, barry (nsf)
c***source             
c***purpose            preconditioner for iterative davidson eigenvalue.
c***                   
c***references         
c
c***routines called    
c***end prologue       trmat
      subroutine trmat(hbuf,ibuf,diag,v,ham,eig,u,t1,trials,etrial,
     1                 scale,dim,n,ntrials,lenbuf,incore)
      implicit integer (a-z)
      real*8 hbuf, diag, v, ham, eig, u, t1, tmp, scale
      real*8 trials, etrial
      character*1 itoc
      character*80 title
      logical incore
      dimension hbuf(lenbuf), ibuf(2,lenbuf), diag(n), v(n)
      dimension ham(n,n), eig(n), u(n,n), t1(*)
      dimension trials(n,ntrials), etrial(ntrials)
      common/io/inp, iout
      write(iout,1) n
      if(incore) then
         call iosys('read integer "'//itoc(dim)//'d number of '//
     1              'elements" from ham',1,nel,0,' ')
         do 10 i=1,n
            ham(i,i) = diag(i)
 10      continue  
         do 20 i=1,nel
            ii=ibuf(1,i)
            jj=ibuf(2,i)
            ham(ii,jj) = hbuf(i)
            ham(jj,ii) = hbuf(i)
 20      continue 
      else
         call rdham(hbuf,ibuf,diag,ham,dim,lenbuf,n,title,.true.)
      endif
      do 30 i=1,n
         ham(i,i) = ham(i,i) - v(i)
 30   continue   
      call tred2(n,n,ham,eig,t1,u)
      call tql2(n,n,eig,t1,u,ierr)
      if(ierr.ne.0) then
          call lnkerr('error from direct diagonalization routine')
      endif
      title='trial eigenvalues'
      call prntfm(title,etrial,ntrials,1,ntrials,1,iout)
      title='trial eigenvectors'
      call prntrm(title,trials,n,ntrials,n,ntrials,iout)
      call vscale(eig,eig,scale,n)
      title='eigenvalues of preconditioner'
      call prntfm(title,eig,n,1,n,1,iout)
      call vscale(eig,eig,1.d0/scale,n)
      title='eigenvectors of preconditioner'
      call prntrm(title,u,n,ntrials,n,ntrials,iout)
      do 40 i=1,ntrials
         etrial(i)=eig(i)
         do 50 j=1,n
            trials(j,i)=u(j,i)
 50      continue   
 40   continue   
      return
 1    format(/,5x,'entering preconditioning routine for matrix n = ',i4)
 2    format(/,5x,'root = ',i3,2x,'energy = ',e15.8)
      end       



