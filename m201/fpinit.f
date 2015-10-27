*deck @(#)fpinit.f	5.1  11/6/94
      subroutine fpinit(nvar,nz,ops,fpcycl,mxcycl,newcyc,inscal,exit,
     $                  tstcrv,chkpt,prnt,convf,toang,pool0,pool1,
     $                  vname,delvar,yold,d1var,d2var,d1vold,xi,h,
     $                  lbl,lalpha,lbeta,intvec,fpvec)
c***begin prologue     fpinit.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)fpinit.f	5.1   11/6/94
c***purpose            
c***description
c
c     initialize fletcher-powell routine for geometry optimization.
c***references
c
c***routines called
c
c***end prologue       fpinit.f
      implicit none
c     --- input variables -----
      integer nvar,nz
      character*(*) ops
      real*8 toang
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer lbl(nz),lalpha(nz),lbeta(nz),intvec(nvar)
      real*8 pool0(nvar),pool1(nvar),fpvec(nvar),delvar(nvar)
      real*8 yold(nvar),d1var(nvar),d2var(nvar),d1vold(nvar)
      real*8 xi(nvar),h(nvar*nvar)
c     --- output variables ---
      integer fpcycl,mxcycl
      logical newcyc,inscal,exit
      logical chkpt,tstcrv,prnt
      real*8 convf
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i,j,intkey,nwords
      logical bond,angle,dump,debug
      logical logkey
      real*8 toau,torad,one,f45,delstr,delbnd
      real*8 pt0001
      real*8 fpkey
c
      parameter (one=1.0d+00,f45=4.5d+01,delstr=1.0d-02,delbnd=1.0d+00)
      parameter (pt0001=1.0d-04)
c
      common/io/inp,iout
c
 1000 format(5x,'symbol:',a16,' in fpinit is neither bond length nor ',
     $          'angle.')
c
c     --- edit local options.
c     --- convergence test on the gradient.
      convf=3.0d-04
      convf=fpkey(ops,'opt=convf',convf,' ')
c
c     --- maximum number of cycles before aborting.
      mxcycl=5
      if(convf.le.pt0001) mxcycl=8
      mxcycl=intkey(ops,'opt=cycles',mxcycl,' ')
c
c     --- checkpoint switch.
      chkpt=logkey(ops,'opt=chkpnt',.true.,' ')
c
c     --- abort test on curvature of stationary point.
      tstcrv=logkey(ops,'opt=tstcrv',.true.,' ')
c
c     --- print switches.
      dump=logkey(ops,'opt=dump',.false.,' ')
      debug=logkey(ops,'opt=debug',.false.,' ')
      prnt=dump.or.debug
c
c
      fpcycl=0
      newcyc=.true.
      inscal=.false.
      exit=.true.
c
c     --- move data read by m101 to the fp local arrays.
      call iosys('read real zvalues from rwf after rewinding',nvar,
     $              pool0,0,' ')
      call vmove(pool1,pool0,nvar)
      nwords=nvar*len(vname(1))
      call iosys('read character "variable names" from rwf',nwords,
     $            0,0,vname)
c
c     --- set the step sizes.
      call iosys('read integer zlbl from rwf',nz,lbl,0,' ')
      call iosys('read integer zlalpha from rwf',nz,lalpha,0,' ')
      call iosys('read integer zlbeta from rwf',nz,lbeta,0,' ')
      call iosys('read integer zintvec from rwf',nvar,intvec,0,' ')
      call iosys('read real zfpvec from rwf',nvar,fpvec,0,' ')
      toau=one/toang
      torad=atan(one)/f45
      do 100 i=1,nvar
         bond=.false.
         angle=.false.
         do 80 j=2,nz
            if(abs(lbl(j)).eq.i) then
               bond=.true.
               goto 90
            else if(j.ge.3) then
               if(abs(lalpha(j)).eq.i.or.abs(lbeta(j)).eq.i) then
                  angle=.true.
                  goto 90
               endif
            endif
   80    continue
c
   90    if(bond) then
            delvar(i)=delstr*toau
            if(intvec(i).ne.0) delvar(i)=fpvec(i)*toau
         else if(angle) then
            delvar(i)=delbnd*torad
            if(intvec(i).ne.0) delvar(i)=fpvec(i)*torad
         else
            write(iout,1000) vname(i)
            call lnkerr(' ')
         endif
  100 continue
c
c     --- initalize derivative and second derivative arrays
      call rzero(yold,nvar)
      call rzero(d1var,nvar)
      call rzero(d2var,nvar)
      call rzero(d1vold,nvar)
      call rzero(xi,nvar)
      call rzero(h,nvar*nvar)
c
c
      return
      end
