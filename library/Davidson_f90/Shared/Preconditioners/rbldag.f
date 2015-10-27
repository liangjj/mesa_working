*deck rbldag.f
c***begin prologue     rbldag
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            construct sub-blocks of the time-dependent 
c                      hamiltonian and perform diagonalizations.
c***                   
c***description        sub-blocks of the full hamiltonian are explicitly 
c***                   constructed and diagonalized.  the eigenvalues
c***                   and transformation matrices are used in the 
c***                   iterative solve as a preconditioner.
c***references         
c
c***routines called    
c***end prologue       rbldag
      subroutine rbldag(ham,eig,hx,hy,hz,vxyz,work,ind,nx,ny,nz,
     1                  n,nblck,dim,prn)
      implicit integer (a-z)
      real*8 ham, eig, hx, hy, hz, vxyz, work
      logical prn
      character*80 title
      character*24 str
      character*2 itoc
      dimension ham(*), eig(*), hx(nx,nx), hy(ny,ny), hz(nz,nz)
      dimension vxyz(n), work(5*n), ind(n,*)
      common/io/inp, iout
      call iosys('create real vectors on ham',n*n,0,0,' ')
      ntrip=n/nblck
      left=n-ntrip*nblck
      if(left.gt.0) then
         ntrip=ntrip+1
      else
         left=nblck
      endif
      if(left.gt.0) then
         write(iout,1) ntrip, nblck, left
      else
         write(iout,2) ntrip, nblck
      endif
      call iosys('write integer "number of blocks" to ham',
     1            1,ntrip,0,' ')
      if(dim.eq.1) then
         call setnd1(ind,nx,n)
      elseif(dim.eq.2) then
         call setnd2(ind,nx,ny,n)
      elseif(dim.eq.3) then
          call setnd3(ind,nx,ny,nz,n)
      endif
      cnt=0
      do 10 trp=1,ntrip
         str='block-'//itoc(trp)
         if(trp.eq.ntrip) then
            msize=left
         else
            msize=nblck
         endif
         call iosys('write integer "block size for '//str//'" to ham',
     1               msize,eig,0,' ') 
         if(dim.eq.1) then
            call rham1(ham,hx,vxyz(cnt+1),ind(cnt+1,1),nx,msize,n,prn)
         elseif(dim.eq.2) then
            call rham2(ham,hx,hy,vxyz(cnt+1),ind(cnt+1,1),nx,ny,
     1                 msize,n,prn)
         elseif(dim.eq.3) then
            call rham3(ham,hx,hy,hz,vxyz(cnt+1),ind(cnt+1,1),nx,ny,nz,
     1                 msize,n,prn)
         endif
         call dsyev('v','l',msize,ham,msize,eig(cnt+1),work,
     1               5*msize,info)
         call iosys('write real "transformation matrix for '
     1               //'" to ham',msize*msize,ham,0,' ')
         call blkout(ham,work,cnt,n,msize)
         title='eigenvalues for '//str
         call prntrm(title,eig(cnt+1),msize,1,msize,1,iout)
         if(prn) then
            title='eigenvectors for '//str
            call prntrm(title,ham,msize,msize,msize,msize,iout)
         endif
         cnt=cnt+msize
 10   continue   
      call iosys('write real "trial eigenvalues" to ham',n,eig,0,' ') 
      return
 1    format(/,5x,'number of blocks in preconditioner = ',i4,/,5x,
     1            'size of biggest block              = ',i4,/,5x,
     2            'size of last block                 = ',i4)
 2    format(/,5x,'number of blocks in preconditioner = ',i4,/,5x,
     1            'size of biggest block              = ',i4)
      end       
