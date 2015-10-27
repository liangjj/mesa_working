*deck @(#)d2einit.f	5.1  11/6/94
      subroutine d2einit(nvar,nvv,maxpt,ops,x,vname,toang,
     $                   xc,natoms)
c***begin prologue     d2einit.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             page, michael(nrl)
c***source             @(#)d2einit.f	5.1   11/6/94
c***purpose            initialize the second derivative processing.
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       d2einit.f
      implicit none
c     --- input variables -----
      integer nvar,nvv,maxpt,natoms
      real*8 toang
c     --- input arrays (unmodified) ---
      character*(*) ops(*)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      character*(*) vname(nvar)
      real*8 x(nvar),xc(3*natoms)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer d2ecycl
      integer inp,iout
      integer j,nwords
      integer d2coord
      logical prnt,chkpt,singpt,cartfx
      real*8 energy,rmax,rmin,rlim,stpsize
      logical logkey,debug
      real*8 fpkey
c
      parameter (debug=.false.)
c
      common/d2einf/energy,rmax,rmin,rlim,d2ecycl,
     $               prnt,chkpt,singpt,stpsize,cartfx
      common/io/inp,iout
c
 1000 format(5x,'differencing cartesian coordinate gradients')
 1010 format(' just read coordinates in d2einit ',/
     $      ,(3f20.10))
c
c     --- pick up some link options.
      singpt=.true.
      singpt=.not.logkey(ops,'force-constants=numerical=two-point',
     $                  .false.,' ')
      stpsize=fpkey(ops,'force-constants=numerical=stepsize',
     $              0.01d+00,' ')
      cartfx=.false.
      if(logkey(ops,'force-constants=numerical=cartesian',.false.,' '))
     $   then
         cartfx=.true.
         d2coord=1
         call iosys('write integer d2e_coord to rwf',1,
     $               d2coord,0,' ')
         write(iout,1000)
      endif
c
c
c     --- move data read by m101 into local arrays.
      d2ecycl=0
      nwords=nvar*len(vname(1))
      call iosys('read character "variable names" from rwf',
     $           nwords,0,0,vname)
      call iosys('read real zvalues from rwf',nvar,x,0,' ')
      call iosys('read real coordinates from rwf',natoms*3,xc,0,' ')
      if(debug) then
         write(iout,1010) (xc(j),j=1,natoms*3)
      endif
c
c
      return
      end
