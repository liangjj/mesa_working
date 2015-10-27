*deck @(#)d2einfo.f	5.1  11/6/94
      subroutine d2einfo(action,nvar,nvv,maxpt,vname,x,f,frcnst,xx,
     $                  ff,xc,fc,xxc,ffc,fs,natoms)
c***begin prologue     d2einfo.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)d2einfo.f	5.1   11/6/94
c***purpose            writes and reads information pertainig to second
c                      derivatives.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       d2einfo.f
      implicit none
c     --- input variables -----
      integer nvar,nvv,maxpt,natoms
      character*(*) action
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
      real*8 x(nvar),f(nvar),frcnst(nvv),xx(nvar,maxpt)
      real*8 ff(nvar,maxpt),fs(maxpt)
      real*8 xc(3*natoms),fc(3*natoms)
      real*8 xxc(natoms*3,maxpt),ffc(natoms*3,maxpt)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer d2ecycl
      integer inp,iout,nwords
      integer lend2e,nat3m,wptoin
      logical prnt,chkpt,singpt,cartfx
      logical debug
      real*8 energy,rmax,rmin,rlim,stpsize
c
      parameter (debug=.false.)
c
      common/d2einf/energy,rmax,rmin,rlim,d2ecycl,
     $               prnt,chkpt,singpt,stpsize,cartfx
      common/io/inp,iout
c
 1000 format(' in d2einfo: d2ecycl =',i5)
c
c     --- either read from or write to the rwf.
      lend2e=5+wptoin(5)
      nat3m=3*natoms*maxpt
      if(action.eq.'write') then
         call iosys('write integer d2e_cycle to rwf',1,d2ecycl,
     $              0,' ')
         call iosys('write integer /d2einf/ to rwf',lend2e,energy,0,' ')
         call iosys('write real "d2e x" to rwf',nvar,x,0,' ')
         call iosys('write real forces to rwf',nvar,f,0,' ')
         call iosys('write real force_constants to rwf',nvv,frcnst,
     $               0,' ')
         call iosys('write real "d2e xx" to rwf',nvar*maxpt,xx,0,' ')
         call iosys('write real "d2e ff" to rwf',nvar*maxpt,ff,0,' ')
         call iosys('write real "d2e xc" to rwf',natoms*3,xc,0,' ')
c        call iosys('write real "cartesian first derivatives" to rwf',
c     $                                             natoms*3,fc,0,' ')
         call iosys('write real "d2e xxc" to rwf',nat3m,xxc,0,' ')
         call iosys('write real "d2e ffc" to rwf',nat3m,ffc,0,' ')
         call iosys('write real "d2e fs" to rwf',maxpt,fs,0,' ')
      else if(action.eq.'read') then
         call iosys('read integer d2e_cycle from rwf',1,
     $              d2ecycl,0,' ')
         call iosys('read integer /d2einf/ from rwf',lend2e,energy,
     $               0,' ')
         call iosys('read real energy from rwf',1,energy,0,' ')
         call iosys('read real "d2e x" from rwf',nvar,x,0,' ')
         call iosys('read real forces from rwf',nvar,f,0,' ')
         call iosys('read real force_constants from rwf',nvv,frcnst,
     $              0,' ')
         call iosys('read real "d2e xx" from rwf',nvar*maxpt,xx,0,' ')
         call iosys('read real "d2e ff" from rwf',nvar*maxpt,ff,0,' ')
         call iosys('read real "d2e xc" from rwf',natoms*3,xc,0,' ')
         call iosys('read real "cartesian first derivatives" from rwf',
     $                                            natoms*3,fc,0,' ')
         call iosys('read real "d2e xxc" from rwf',nat3m,xxc,0,' ')
         call iosys('read real "d2e ffc" from rwf',nat3m,ffc,0,' ')
         call iosys('read real "d2e fs" from rwf',maxpt,fs,0,' ')
         nwords=nvar*len(vname(1))
         call iosys('read character "variable names" from rwf',nwords,
     $              0,0,vname)
      endif
c
      if(debug) then
         write(iout,1000) d2ecycl
      endif
c
c
      return
      end
