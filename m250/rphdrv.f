*deck @(#)rphdrv.f	5.1  11/6/94
      subroutine rphdrv(nvar,nz,natoms,ops,vname,ian,atmass,
     $                  cmass,xc,gc,zf,ftric,fsq,nvvc,order,
     $                  rphdone,rphstep,stotal)
c***begin prologue     rphdrv.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)rphdrv.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       rphdrv.f
      implicit none
c     --- input variables -----
      integer nvar,nz,natoms,nvvc
c     --- input arrays (unmodified) ---
      integer ian(nz),order(3*natoms)
      character*(*) ops,vname(nvar)
      real*8 atmass(natoms)
      real*8 ftric(nvvc),cmass(3*natoms)
c     --- input arrays (scratch) ---
      real*8 gc(3*natoms)
      real*8 fsq(3*natoms,3*natoms)
      real*8 zf(*),xc(3*natoms)
c     --- output arrays ---
c     --- output variables ---
      integer rphstep
      logical rphdone
      real*8 sstep
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer intkey,nat3,i,rphmax,iphase
      character*4 ians
      logical projg
      real*8 stotal,gnorm
      real*8 fpkey,phase,zero,one
      logical logkey,saddle
c
      parameter (zero=0.0d+00,one=1.0d+00)
c
      common/io/inp,iout
c
 1000 format(/' m250:        step number = ',i4,
     $       /'                 stepsize = ',f6.3,
     $       /'               arc length = ',f6.3,
     $       /'  maximum number of steps = ',i4)
 1010 format(/'     norm of the gradient = ',f10.6)
 1020 format(/' beginning a steepest descent reaction pathway',
     $       /, ' at a non stationary point')
 1030 format(/' beginning a steepest descent reaction pathway',
     $       /, ' at a stationary point')
c
c     --- see if this is the first step
      call iosys('does rph_step exist on rwf',0,0,0,ians)
      if(ians.eq.'no') then
         rphstep=0
         stotal=zero
      else
         call iosys('read integer rph_step from rwf',1,rphstep,0,' ')
         call iosys('read real rph_stotal from rwf',1,stotal,0,' ')
      endif
c
      if(rphstep.eq.1) then
         call iosys('write character "print flag" to rwf',
     $               0,0,0,'minimum')
      endif
c
      rphmax=intkey(ops,'rph=maximum-number-of-steps',5,' ')
      sstep=fpkey(ops,'rph=stepsize',0.01d+00,' ')
      iphase=intkey(ops,'rph=phase',1,' ')
      saddle=logkey(ops,'rph=saddle',.false.,' ')
      if(iphase.eq.1) then
         phase=one
      else if(iphase.eq.-1) then
         phase=-one
      else
         call lnkerr(' m250: dont understand phase, try =1 or =-1')
      endif
c
      write(iout,1000) rphstep,sstep,stotal,rphmax
c
      if(rphstep.ne.0) then
         nat3=natoms*3
         call iosys('read real "cartesian first derivatives" from rwf',
     $              nat3,gc,0,' ')
         gnorm=zero
         do 100 i=1,nat3
            gnorm=gnorm+gc(i)**2
  100    continue
         gnorm=sqrt(gnorm)
         write(iout,1010) gnorm
         if(stotal.eq.0.0) then
c
c           --- first real entry
c               is this a stationary point?
            if(gnorm.gt.1.0d-3) then
               write(iout,1020)
               projg=.true.
            else
               write(iout,1030)
               projg=.false.
            endif
         else
c
c           --- subsequent entry
            projg=.true.
         endif
c
c
         call frmain(nvar,nz,natoms,ops,vname,ian,atmass,
     $               cmass,xc,zf,ftric,fsq,nvvc,order,projg,
     $               sstep,phase,saddle)
c
         stotal=stotal+sstep
      endif
c
c     --- increment the step counter
      rphstep=rphstep+1
      if(rphstep.ge.rphmax) rphdone=.true.
c
c
      return
      end
