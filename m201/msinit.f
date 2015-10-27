*deck @(#)msinit.f	5.1  11/6/94
      subroutine msinit(nvar,nvv,ops,zsq,fc,frcnst,id2e,x,vname,hinv,
     $                  nz,lbl,toang,srcd2e,ianz,iz,lalpha,lbeta,
     $                  xangst,isave)
c***begin prologue     msinit.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkey, et al. (g82)
c***source             @(#)msinit.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       msinit.f
      implicit none
c     --- input variables -----
      integer nvar,nvv,nz
c     --- input arrays (unmodified) ---
      character*(*) ops,vname(nvar),srcd2e
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer id2e(nvar),lbl(nz)
      integer ianz(nz),iz(4,nz),lalpha(nz),lbeta(nz),isave(nvar)
      real*8 zsq(nvar,nvar),fc(nvv),x(nvar),frcnst(nvv),hinv(nvar,nvar)
      real*8 xangst(nvar)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer mxcycl,mscycl,np,neg,ipsav1,ipsav2
      integer icalcf,icalff,irdfcd
      integer i,j,ij,idx,intkey,lind,nwords
      character*8 chrkey
      logical prnt,tstcrv,dump,debug,logkey
      logical updrwf,chkpt,clnup
      real*8 convf,dxmaxt,fmaxt,eigmax,eigmin,delta,cutoff
      real*8 energy,rmax,rmin,rlim,fncerr,grderr,fnccnv,alpha
      real*8 det,zero,ten,toang,pt0001
      real*8 fpkey
c
      parameter (zero=0.0d+00,ten=1.0d+01)
      parameter (icalcf=3,icalff=4,irdfcd=1)
      parameter (pt0001=1.0d-04)
      parameter (cutoff=1.0d-06,delta=5.0d-03)
c
      common/io/inp,iout
      common/msinf/ energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,
     $              eigmax,eigmin,fncerr,grderr,fnccnv,alpha,
     $              mxcycl,mscycl,np,neg,ipsav1,ipsav2,
     $              prnt,tstcrv,updrwf,chkpt,clnup
c
 1000 format(5x,'expert switch is set.  fmaxt,dxmaxt,eigmax,eigmin',
     $          'relaxed.')
 1010 format(8f10.6)
c
      lind(i,j)=(0.5*max(i,j)*(max(i,j)-1)) +min(i,j)
c
c     --- set up some optimization defaults ---
      fncerr=1.0d-07
      grderr=1.0d-07
      fnccnv=1.0d+00
      rmax=6.0d-01
      rmin=1.5d-03
      rlim=7.0d-02
c
c     --- edit local options list.
c         --- curvature requirement for the stationary point. 
c             default to local minimum.
      if(logkey(ops,'opt=ts',.false.,' ')) then
         neg=1
      else
         neg=intkey(ops,'opt=negeig',0,' ')
      endif
c
c     --- convergence on the gradient.
      convf=3.0d-04
      convf=fpkey(ops,'opt=convf',convf,' ')
c
c     --- maximum number of cycles before aborting.
      mxcycl=min(20,nvar+10)
      if(convf.le.pt0001) mxcycl=3*mxcycl/2
      mxcycl=intkey(ops,'opt=cycles',mxcycl,' ')
c
c     --- maximum step size allowed.
      dxmaxt=2.0d-01
      dxmaxt=fpkey(ops,'opt=dxmaxt',dxmaxt,' ')
c
c     --- source of initial second derivatives.
      srcd2e=chrkey(ops,'opt=srcd2e',' ',' ')
c
c     --- abort test on curvature of stationary point.
      tstcrv=logkey(ops,'opt=tstcrv',.true.,' ')
c
c     --- abort test for big forces.
      fmaxt=fpkey(ops,'opt=fmaxt',1.0d+00,' ')
c
c     --- acceptable eigenvalue ranges for second derivative matrix.
      eigmax=2.5d+01
      eigmin=1.0d-04
      eigmax=fpkey(ops,'opt=eigmax',eigmax,' ')
      eigmin=fpkey(ops,'opt=eigmin',eigmin,' ')
c
c     --- expert switch.
      if(logkey(ops,'opt=expert',.false.,' ')) then
         write(iout,1000)
         fmaxt=fmaxt*ten
         dxmaxt=dxmaxt*ten
         eigmax=eigmax*10
         eigmin=eigmin/ten
      endif
c
c     --- rwf update switch.
      updrwf=.true.
c
c     --- print switches.
      dump=logkey(ops,'opt=dump',.false.,' ')
      debug=logkey(ops,'opt=debug',.false.,' ')
      prnt=dump.or.debug
c
c     --- checkpoint switch.
      chkpt=.true.
      chkpt=logkey(ops,'opt=chkpnt',.true.,' ')
c
c     --- d2e file cleanup switch.
      clnup=logkey(ops,'opt=clnd2e',.false.,' ')
c
c     --- set optimization parameters.
      np=0
      mscycl=0
      ipsav1=0
      ipsav2=0
c
c     --- move data read by m101 into local arrays.
      nwords=nvar*len(vname(1))
      call iosys('read character "variable names" from rwf',nwords,
     $            0,0,vname)
      call iosys('read integer zintvec from rwf',nvar,id2e,0,' ')
c     --- if diagonal second derivatives were read by m101, they appear
c         in fpvec.
      call iosys('read real zfpvec from rwf',nvar,x,0,' ')
      call rzero(zsq,nvar*nvar)
      do 10 i=1,nvar
         zsq(i,i)=x(i)
   10 continue
      call iosys('read real zvalues from rwf',nvar,x,0,' ')
c
c     --- initialize the second derivative matrix.
c         the array id2e will tell estmd2 what to do.
c         if id2e(i)=0,  estmd2 should provide an initial guess.
c         if id2e(i)=1,  diagonal initial guess was read by m101.
c         if id2e(i)=2,  full second derivative matrix read by m101
c         if id2e(i)=3,  the derivatives will be estimated numerically.
c         if id2e(i)=4,  the second derivatives will be calculated analytically.
      call rzero(fc,nvv)
      if(srcd2e.eq.' ') then
         do 20 i=1,nvar
            ij=i*(i+1)/2
            fc(ij)=zsq(i,i)
   20    continue
      else if(srcd2e.eq.'analytic') then
         do 30 i=1,nvar
            id2e(i)=icalff
   30    continue
         call rzero(zsq,nvar*nvar)
      else if(srcd2e.eq.'chk'.or.srcd2e.eq.'rdchk') then
         write(iout,*) ' initial force constants obtained from chk file'
         call readfc(frcnst,nvv,'chk')
         call vmove(fc,frcnst,nvv)
         do 40 i=1,nvar
            idx=i*(i+1)/2
            if(id2e(i).eq.0) then
               id2e(i)=1
            else if(id2e(i).eq.1) then
               fc(idx)=zsq(i,i)
            else
               if(zsq(i,i).ge.cutoff) then
                  fc(idx)=zsq(i,i)
               else
                  fc(idx)=delta
               endif
            endif
   40    continue
      else if(srcd2e.eq.'rdinp') then
c        read(inp,1010) ((zsq(i,j),j=1,i),i=1,nvar)
c        do 50 i=1,nvar
c           if(id2e(i).eq.icalcf) then
c              if(abs(zsq(i,i)).le.cutoff) zsq(i,i)=delta
c           else
c              if(abs(zsq(i,i)).ge.cutoff) id2e(i)=irdfcd
c           endif
c50      continue
         write(iout,*) 'initial force constants obtained from input'
         call readfc(frcnst,nvv,'rwf')
         call vmove(fc,frcnst,nvv)
         do 50 i=1,nvar
            idx=i*(i+1)/2
            if(id2e(i).eq.0) then
               id2e(i)=2
            else if(id2e(i).eq.1) then
               fc(idx)=zsq(i,i)
            else
               if(zsq(i,i).ge.cutoff) then
                  fc(idx)=zsq(i,i)
               else
                  fc(idx)=delta
               endif
            endif
 50      continue
      endif
c
c     --- obtain preliminary set of diagonal second derivatives.
c         this section is processed even if second derivatives were read
c         from input.
      if(srcd2e.ne.'analytic') then
c        --- retrieve the z-matrix information.
         call iosys('read integer zian from rwf',nz,ianz,0,' ')
         call iosys('read integer ziz from rwf',4*nz,iz,0,' ')
         call iosys('read integer zlalpha from rwf',nz,lalpha,0,' ')
         call iosys('read integer zlbeta from rwf',nz,lbeta,0,' ')
c
         call estmd2(nvar,nvv,nz,toang,x,fc,id2e,ianz,iz,lbl,lalpha,
     $             lbeta,xangst,isave)
         do 60 i=1,nvar
            do 60 j=1,i
               if((id2e(i).eq.icalff).and.(id2e(j).eq.icalff))
     $                     fc(lind(i,j))=frcnst(lind(i,j))
   60    continue
      endif
c
c     --- print information concerning initial second derivatives.
      call prmtbl(0,vname,x,id2e,fc,nvar,lbl,nz,toang)
c
c     --- form initial approximation to the inverse hessian.
c         if analytic derivatives are to be used,this is unnecessary.
      if(srcd2e.ne.'analytic') then
         call trtosq(hinv,fc,nvar,nvv)
         call minvrt(hinv,nvar,nvar,det,isave,xangst)
      else
         call rzero(hinv,nvar*nvar)
      endif
c
c
      return
      end
