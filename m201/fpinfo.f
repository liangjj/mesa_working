*deck @(#)fpinfo.f	5.1  11/6/94
      subroutine fpinfo(action,nvar,pool0,pool1,delvar,yold,
     $                  d1var,d2var,d1vold,xi,h,vname)
c***begin prologue     fpinfo.f
c***date written       yymmdd  
c***revision date      11/6/94      
c   18 october 1992, rlm at lanl
c      fixing bug with lenfp on 32/64 bit machine
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)fpinfo.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       fpinfo.f
      implicit none
c     --- input variables -----
      integer nvar
      character*(*) action
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
      real*8 pool0(nvar),pool1(nvar),delvar(nvar),yold(nvar)
      real*8 d1var(nvar),d2var(nvar),d1vold(nvar),xi(nvar),h(nvar,nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer fpcycl,mxcycl,mode,nstep,curvar,lambda
      integer lenfp,nwords,wptoin
      real*8 convf,alpha,fzero,f1
      logical newcyc,inscal
c
      common/io/inp,iout
      common/fpinf/ fpcycl,mxcycl,convf,mode,nstep,curvar,lambda,newcyc,
     $              inscal,alpha,fzero,f1(4)
c
c     --- determine length of common block
      lenfp=wptoin(7)+8
      if(action.eq.'write') then
         call iosys('write integer "optimization cycle" to rwf',
     $               1,fpcycl,0,' ')
         call iosys('write integer /fpinfo/ to rwf',lenfp,fpcycl,0,' ')
         call iosys('write real "fp pool0" to rwf' ,nvar,pool0,0,' ')
         call iosys('write real "fp pool1" to rwf' ,nvar,pool1,0,' ')
         call iosys('write real "fp delvar" to rwf' ,nvar,delvar,0,' ')
         call iosys('write real "fp yold" to rwf' ,nvar,yold,0,' ')
         call iosys('write real "fp d1var" to rwf' ,nvar,d1var,0,' ')
         call iosys('write real "fp d2var" to rwf',nvar,d2var,0,' ')
         call iosys('write real "fp d1vold" to rwf',nvar,d1vold,0,' ')
         call iosys('write real "fp xi" to rwf',nvar,xi,0,' ')
         call iosys('write real "fp h" to rwf',nvar*nvar,h,0,' ')
         nwords=nvar*len(vname(1))
         call iosys('write character "variable names" to rwf',nwords,0,
     $              0,vname)
c        --- update the z-matrix substitution values.
         call iosys('write real zvalues to rwf',nvar,pool0,0,' ')
      else if(action.eq.'read') then
         call iosys('read integer /fpinfo/ from rwf',lenfp,fpcycl,0,' ')
         call iosys('read real "fp pool0" from rwf' ,nvar,pool0,0,' ')
         call iosys('read real "fp pool1" from rwf' ,nvar,pool1,0,' ')
         call iosys('read real "fp delvar" from rwf' ,nvar,delvar,0,' ')
         call iosys('read real "fp yold" from rwf' ,nvar,yold,0,' ')
         call iosys('read real "fp d1var" from rwf' ,nvar,d1var,0,' ')
         call iosys('read real "fp d2var" from rwf',nvar,d2var,0,' ')
         call iosys('read real "fp d1vold" from rwf',nvar,d1vold,0,' ')
         call iosys('read real "fp xi" from rwf',nvar,xi,0,' ')
         call iosys('read real "fp h" from rwf',nvar*nvar,h,0,' ')
         nwords=nvar*len(vname(1))
         call iosys('read character "variable names" from rwf',nwords,
     $              0,0,vname)
      endif
c
c
      return
      end
