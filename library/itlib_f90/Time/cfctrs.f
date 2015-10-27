*deck cfctrs.f
c***begin prologue     cfctrs
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            construct sub-blocks of the time-dependent 
c                      hamiltonian and perform LU factorizations.
c***                   
c***description        sub-blocks of the full time-dependent
c***                   hamiltonian in complex form are explicitly 
c***                   constructed and the LU factorizations are performed.
c***                   the factorization is used in the iterative solve 
c***                   as a preconditioner.  the blocks are not allowed to
c***                   cross channel boundaries by construction.
c***references         
c
c***routines called    
c***end prologue       cfctrs
      subroutine cfctrs(ham,driver,vtrial,tmp,hx,hy,hz,ht,vxyzt,
     1                   ind,ilst,nx,ny,nz,nt,nc,n,nblck,dim,
     2                   prn,trial)
      implicit integer (a-z)
      real*8 hx, hy, hz, vxyzt, driver
      real*8 ht, vtrial
      complex*16 ham, tmp
      logical prn, trial
      character*80 title
      character*2 itoc
      dimension ham(*), ht(nt,nt), vxyzt(n,nc,nc)
      dimension ind(n,*), ilst(nblck)
      dimension hx(nx,nx), hy(ny,ny), hz(nz,nz)
      dimension driver(n,nc,2), tmp(n,nc), vtrial(n,nc,2)
      common/io/inp, iout      
c
c     calculate the number of passes needed to block factorize
c     the sub-matrix of the hamiltonian consisting of one super-block
c     which is diagonal in channels.
c
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
c
c      the total number of blocks is nc times ntrip
c
      call iosys('write integer "number of LU blocks '//
     1           'per channel" to ham',1,ntrip,0,' ')
c
c     set up the global index array to label the matrix elements for
c     the full hamiltonian
c
      if(dim.eq.1) then
         call setnd2(ind,nx,nt,n)
      elseif(dim.eq.2) then
         call setnd3(ind,nx,ny,nt,n)
      elseif(dim.eq.3) then
          call setnd4(ind,nx,ny,nz,nt,n)
      endif
c
c     convert the rhs to a complex vector
c
      if(trial) then
         call rv2cv(driver,tmp,n,nc)
      endif
c
c     perform the factorization channel by channel
c
      do 10 ic=1,nc 
         cnt=0
         do 20 i=1,ntrip
            if(i.eq.ntrip) then
               msize=left
            else
               msize=nblck
            endif
c
c           in the calls to cham(1,2,3) the potential is treated
c           as being totally diagonal.  since we are not crossing channel
c           boundaries, that is the case.
c  
            if(dim.eq.1) then
               call cham1(ham,hx,ht,vxyzt(cnt+1,ic,ic),
     1                    ind(cnt+1,1),nx,nt,msize,n,prn)
            elseif(dim.eq.2) then
               call cham2(ham,hx,hy,ht,vxyzt(cnt+1,ic,ic),
     1                    ind(cnt+1,1),nx,ny,nt,msize,n,prn)
            elseif(dim.eq.3) then
               call cham3(ham,hx,hy,hz,ht,vxyzt(cnt+1,ic,ic),
     1                    ind(cnt+1,1),nx,ny,nz,nt,msize,n,prn)
            endif
c
c
c           factorize the complex sub-block and perform the solve.
c
            call csublu(ham,tmp(cnt+1,ic),ilst,msize,prn,trial)
c
c
c           write out the information to the disk.
c
            call iosys('write integer "block size for block-'//itoc(i)
     1                 //'channel-'//itoc(ic)//'" to ham',1,msize,0,' ')
            call iosys('write real "LU factors for block-'//itoc(i)
     1                  //'channel-'//itoc(ic)//'" to ham',
     2                  2*msize*msize,ham,0,' ')
            call iosys('write integer "LU array for block-'//itoc(i)
     1                 //'channel-'//itoc(ic)//'" to ham',
     2                 msize,ilst,0,' ') 
            cnt=cnt+msize
 20      continue   
 10   continue   
      if(trial) then
         call cv2rv(vtrial,tmp,n,nc)
         if(prn) then
            title='trial solution'
            call prntcm(title,tmp,n*nc,1,n*nc,1,iout)
         endif
         call iosys('write real "LU trial vector" to ham',
     1               n*nc*2,vtrial,0,' ') 
      endif
      return
 1    format(/,5x,'number of blocks in preconditioner = ',i4,/,5x,
     1            'size of biggest block              = ',i4,/,5x,
     2            'size of last block                 = ',i4)
 2    format(/,5x,'number of blocks in preconditioner = ',i4,/,5x,
     1            'size of biggest block              = ',i4)
      end       
