*deck rfctrs.f
c***begin prologue     rfctrs
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           time development
c***author             schneider, barry (nsf)
c***source             
c***purpose            construct sub-blocks of the time-dependent 
c                      hamiltonian and perform LU factorizations.
c***                   
c***description        explicit sub-blocks of the full time-dependent
c***                   hamiltonian in complex form are explicitly 
c***                   constructed and the LU factorizations are performed.
c***                   the factorization is used in the iterative solve 
c***                   as a preconditioner.
c***references         
c
c***routines called    
c***end prologue       rfctrs
      subroutine rfctrs(ham,driver,vtrial,tmp,hx,hy,hz,vxyz,
     1                   ind,ilst,nx,ny,nz,nc,n,nblck,dim,
     2                   prn,trial)
      implicit integer (a-z)
      real*8 hx, hy, hz, vxyz, driver
      real*8 vtrial, ham, tmp
      logical prn, trial
      character*80 title
      character*2 itoc
      dimension ham(*), vxyz(n,nc,nc), ind(n,*), ilst(nblck)
      dimension hx(nx,nx), hy(ny,ny), hz(nz,nz)
      dimension driver(n,nc), tmp(n,nc), vtrial(n,nc)
      common/io/inp, iout      
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
      call iosys('write integer "number of LU blocks '//
     1           'per channel" to ham',1,ntrip,0,' ')
      if(dim.eq.1) then
         call setnd1(ind,nx,n)
      elseif(dim.eq.2) then
         call setnd2(ind,nx,ny,n)
      elseif(dim.eq.3) then
          call setnd3(ind,nx,ny,nz,n)
      endif
      do 10 ic=1,nc
         do 20 i=1,ntrip
            if(i.eq.ntrip) then
               msize=left
            else
               msize=nblck
            endif
            if(dim.eq.1) then
               call rham1(ham,hx,vxyz(cnt+1,ic,ic),
     1                    ind(cnt+1,1),nx,msize,n,prn)
            elseif(dim.eq.1) then
               call rham2(ham,hx,hy,vxyz(cnt+1,ic,ic),
     1                    ind(cnt+1,1),nx,ny,msize,n,prn)
            elseif(dim.eq.3) then
               call rham3(ham,hx,hy,hz,vxyz(cnt+1,ic,ic),
     1                    ind(cnt+1,1),nx,ny,nz,msize,n,prn)
            endif
            call rsublu(ham,tmp(cnt+1,ic),ilst,msize,prn,trial)
            call iosys('write integer "block size for block-'//itoc(i)
     1                 //'channel-'//itoc(ic)//'" to ham',1,msize,0,' ')
            call iosys('write real "LU factors for block-'//itoc(i)
     1                  //'channel-'//itoc(ic)//'" to ham',
     2                  msize*msize,ham,0,' ')
            call iosys('write integer "LU array for block-'//itoc(i)
     1                 //'channel-'//itoc(ic)//'" to ham',
     2                 msize,ilst,0,' ') 
            cnt=cnt+msize
 20      continue   
 10   continue   
      if(trial) then
         if(prn) then
            title='solution'
            call prntrm(title,vtrial,n*nc,1,n*nc,1,iout)
         endif
         call iosys('write real "LU trial vector" to ham',
     1               n*nc,vtrial,0,' ') 
      endif
      return
 1    format(/,5x,'number of blocks in preconditioner = ',i4,/,5x,
     1            'size of biggest block              = ',i4,/,5x,
     2            'size of last block                 = ',i4)
 2    format(/,5x,'number of blocks in preconditioner = ',i4,/,5x,
     1            'size of biggest block              = ',i4)
      end       
