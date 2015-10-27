*deck @(#)fpmain.f	5.2  4/18/95
      subroutine fpmain(ops,vname,nvar,nz,convrg,abnrml,toang,pool0,
     $                  pool1,delvar,yold,d1var,d2var,d1vold,xi,h,
     $                  lbl,lalpha,lbeta,intvec,fpvec,z)
c***begin prologue     fpmain.f
c***date written       yymmdd  
c***revision date      4/18/95      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)fpmain.f	5.2   4/18/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       fpmain.f
      implicit none
c     --- input variables -----
      integer nvar,nz
      character*(*) ops
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
c     --- input arrays (scratch) ---
      integer lbl(nz),lalpha(nz),lbeta(nz)
      integer intvec(nvar)
      real*8 pool0(nvar),pool1(nvar),delvar(nvar),yold(nvar)
      real*8 d1var(nvar),d2var(nvar),d1vold(nvar),xi(nvar)
      real*8 h(nvar,nvar),fpvec(nvar),z(*)
c     --- output arrays ---
c     --- output variables ---
      real*8 toang
      logical convrg,abnrml
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer fpcycl,mxcycl,mode,nstep,curvar,lambda
      integer i,intkey
      character*128 namchk
      character*4 ians
      real*8 convf,alpha,fzero,f1,f
      logical dump,exit,inscal,newcyc,archiv,chkpt
      logical logkey
      logical tstcrv
c
      data dump/.false./
      data chkpt/.true./
      save dump,chkpt
c
      common/fpinf/ fpcycl,mxcycl,convf,mode,nstep,curvar,lambda,newcyc,
     $              inscal,alpha,fzero,f1(4)
      common/io/inp,iout
c
 1000 format(1x,'fletcher-powell optimization:')
 1010 format(5x,'initialization pass.')
 1020 format(5x,'restart from: ',a8)
 1030 format(5x,a16,f9.6)
 1040 format(5x,'energy:',f18.9)
c
c
      write(iout,1000)
c
c     --- determine if this is a restart, the initial entry of a new
c         optimization, or a continuation of an optimization.  unless this
c         is a continuation, the cycle counter will be missing from the rwf.
      call iosys('does "optimization cycle" exist on rwf',0,0,0,ians)
      if (ians.eq.'no') then
c        --- either first point or a restart.
         if(logkey(ops,'opt=restart',.false.,' ')) then
c           --- initial entry, restart optimiztion.
            call iosys('read character "checkpoint filename" from rwf',
     $                 -1,0,0,namchk)
            write(iout,1020) namchk
            call chkpnt('restore')
            call fpinfo('read',nvar,pool0,pool1,delvar,yold,d1var,d2var,
     $                  d1vold,xi,h,vname)
            mxcycl=intkey(ops,'opt=cycles',mxcycl,' ')
            if(dump) call fpdump(nvar,pool0,pool1,delvar,yold,d1var,
     $                           d2var,d1vold,xi,h,vname)
         else
c           --- initial entry, new optimization.
            write(iout,1010)
            call fpinit(nvar,nz,ops,fpcycl,mxcycl,newcyc,inscal,exit,
     $                  tstcrv,chkpt,dump,convf,toang,pool0,pool1,
     $                  vname,delvar,yold,d1var,d2var,d1vold,xi,h,
     $                  lbl,lalpha,lbeta,intvec,fpvec)
            call fpinfo('write',nvar,pool0,pool1,delvar,yold,d1var,
     $                  d2var,d1vold,xi,h,vname)
            return
         endif
      else
c        --- a continuation of the current run.
         call iosys('read integer "optimization cycle" from rwf',
     $              1,fpcycl,0,' ')
         call iosys('write character "print flag" to rwf',0,0,0,
     $              'minimum')
         call iosys('read real energy from rwf',1,f,0,' ')
         call fpinfo('read',nvar,pool0,pool1,delvar,yold,d1var,d2var,
     $               d1vold,xi,h,vname)
c        --- print current variables and the energy
         do 30 i=1,nvar
            write(iout,1030) vname(i),pool0(i)
   30    continue
         write(iout,1040) f
c
         if(dump) call fpdump(nvar,pool0,pool1,delvar,yold,d1var,
     $                        d2var,d1vold,xi,h,vname)
         if(chkpt) call chkpnt('save')
      endif
c
c     --- let's go.
      if(inscal) then
         call fpsclt(nvar,pool0,pool1,d2var,xi,f,fzero,f1,lambda,
     $               exit)
         if(exit) then
            newcyc=.true.
            inscal=.false.
         endif
      else
         call fpderv(nvar,newcyc,exit,fpcycl,curvar,nstep,mode,dump,
     $               tstcrv,pool0,pool1,d1var,d2var,yold,d1vold,
     $               delvar,xi,h,vname,fzero,f,f1)
         newcyc=.false.
         if(exit) then
c
c           --- have completed a cycle.  are we at the minimum?
            call fptest(nvar,archiv,convrg,abnrml,fpcycl,mxcycl,convf,
     $                  d1var,pool1,yold)
            if(convrg) then
               call fpexit(nvar,nz,toang,lbl,pool0,delvar,z,z(nvar+1),
     $                     vname,fzero)
            else
               call fphmat(nvar,fpcycl,dump,vname,pool0,pool1,d1var,
     $                     d2var,xi,h,yold,d1vold,z,z(nvar+1),
     $                     z(2*nvar+1),z(2*nvar+1+nvar*nvar))
               fpcycl=fpcycl+1
               lambda=1
               inscal=.true.
            endif
         endif
      endif
c
c     --- save the results and return.
      call fpinfo('write',nvar,pool0,pool1,delvar,yold,d1var,d2var,
     $            d1vold,xi,h,vname)
c
c
      return
      end
