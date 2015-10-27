*deck @(#)msinfo.f	5.1  11/6/94
      subroutine msinfo(action,nvar,nvv,mxpt,vname,x,f,fc,frcnst,xx,
     $                  ff,fs,ic,srcd2e,hinv)
c***begin prologue     msinfo.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)msinfo.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       msinfo.f
      implicit none
c     --- input variables -----
      integer nvar,nvv,mxpt
      character*(*) action
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer ic(mxpt)
      real*8 x(nvar),f(nvar),fc(2*nvv),frcnst(nvv),xx(nvar,mxpt)
      real*8 ff(nvar,mxpt),fs(mxpt),hinv(nvar,nvar)
c     --- output variables ---
      character*(*) srcd2e
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer mxcycl,mscycl,np,neg,ipsav1,ipsav2
      integer lenms,nwords,wptoin
      real*8 energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,eigmax,eigmin
      real*8 fncerr,grderr,fnccnv,alpha
      logical prnt,tstcrv,updrwf,chkpt,clnup
c
      common/msinf/ energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,
     $              eigmax,eigmin,fncerr,grderr,fnccnv,alpha,
     $              mxcycl,mscycl,np,neg,ipsav1,ipsav2,
     $              prnt,tstcrv,updrwf,chkpt,clnup
      common/io/inp,iout
c
c     --- determine how long the common block is
      lenms=wptoin(13)+11
      if(action.eq.'write') then
         call iosys('write integer "optimization cycle" to rwf',1,
     $              mscycl,0,' ')
         call iosys('write character source_d2e to rwf',0,0,0,srcd2e)
         call iosys('write integer /msinfo/ to rwf',lenms,energy,0,' ')
         call iosys('write real "ms x" to rwf',nvar,x,0,' ')
         call iosys('write real forces to rwf',nvar,f,0,' ')
         call iosys('write real "ms fc" to rwf',2*nvv,fc,0,' ')
         call iosys('write real force_constants to rwf',nvv,frcnst,
     $               0,' ')
         call iosys('write real "ms xx" to rwf',nvar*mxpt,xx,0,' ')
         call iosys('write real "ms ff" to rwf',nvar*mxpt,ff,0,' ')
         call iosys('write real "ms fs" to rwf',mxpt,fs,0,' ')
         call iosys('write integer "ms ic" to rwf',mxpt,ic,0,' ')
         call iosys('write real "ms inv hessian" to rwf',nvar*nvar,
     $              hinv,0,' ')
         nwords=nvar*len(vname(1))
         call iosys('write character "variable names" to rwf',nwords,
     $              0,0,vname)
      else if(action.eq.'read') then
         call iosys('read integer "optimization cycle" from rwf',1,
     $              mscycl,0,' ')
         call iosys('read character source_d2e from rwf',-1,0,0,srcd2e)
         call iosys('read integer /msinfo/ from rwf',lenms,energy,0,' ')
         call iosys('read real energy from rwf',1,energy,0,' ')
         call iosys('read real "ms x" from rwf',nvar,x,0,' ')
         call iosys('read real forces from rwf',nvar,f,0,' ')
         call iosys('read real "ms fc" from rwf',2*nvv,fc,0,' ')
         call iosys('read real force_constants from rwf',nvv,frcnst,
     $              0,' ')
         call iosys('read real "ms xx" from rwf',nvar*mxpt,xx,0,' ')
         call iosys('read real "ms ff" from rwf',nvar*mxpt,ff,0,' ')
         call iosys('read real "ms fs" from rwf',mxpt,fs,0,' ')
         call iosys('read integer "ms ic" from rwf',mxpt,ic,0,' ')
         call iosys('read real "ms inv hessian" from rwf',nvar*nvar,
     $              hinv,0,' ')
         nwords=nvar*len(vname(1))
         call iosys('read character "variable names" from rwf',nwords,
     $              0,0,vname)
      endif
c
c
      return
      end
