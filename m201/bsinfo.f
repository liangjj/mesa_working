*deck @(#)bsinfo.f	5.1  11/6/94
      subroutine bsinfo(action,nvar,nvv,mxpt,vname,x,f,fc,frcnst,xx,
     $                  ff,fs,ic,srcd2e,srcht)
c***begin prologue     bsinfo.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)bsinfo.f	5.1   11/6/94
c***purpose            reads or writes informations regarding berny optimizations
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       bsinfo.f
      implicit none
c     --- input variables -----
      integer nvar,nvv,mxpt
      character*(*) action
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      integer ic(mxpt)
      character*(*) srcd2e,srcht
      real*8 x(nvar),f(nvar),fc(2*nvv),frcnst(nvv),xx(nvar,mxpt)
      real*8 ff(nvar,mxpt),fs(mxpt)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer mxcycl,bscycl,np,neg,ipsav1,ipsav2
      integer lenbs,nwords,wptoin
      real*8 energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,eigmax,eigmin
      real*8 fswtch,fncerr,grderr,fnccnv
      logical prnt,tstcrv,updrwf,chkpt,clnup
c
      common/bsinf/ energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,
     $              eigmax,eigmin,fswtch,fncerr,grderr,fnccnv,
     $              mxcycl,bscycl,np,neg,ipsav1,ipsav2,
     $              prnt,tstcrv,updrwf,chkpt,clnup
      common/io/inp,iout
c
c     --- get length of common block /bsinf/
      lenbs=wptoin(13)+11
      if(action.eq.'write') then
         call iosys('write integer "optimization cycle" to rwf',
     $        1,bscycl,0,' ')
         call iosys('write character search_type to rwf',0,0,0,srcht)
         call iosys('write character source_d2e to rwf',0,0,0,srcd2e)
         call iosys('write integer /bsinfo/ to rwf',lenbs,energy,0,' ')
         call iosys('write real "bs x" to rwf',nvar,x,0,' ')
         call iosys('write real forces to rwf',nvar,f,0,' ')
         call iosys('write real "bs fc" to rwf',2*nvv,fc,0,' ')
         call iosys('write real force_constants to rwf',nvv,frcnst,
     $               0,' ')
         call iosys('write real "bs xx" to rwf',nvar*mxpt,xx,0,' ')
         call iosys('write real "bs ff" to rwf',nvar*mxpt,ff,0,' ')
         call iosys('write real "bs fs" to rwf',mxpt,fs,0,' ')
         call iosys('write integer "bs ic" to rwf',mxpt,ic,0,' ')
         nwords=nvar*len(vname(1))
         call iosys('write character "variable names" to rwf',nwords,0,
     $              0,vname)
      else if(action.eq.'read') then
         call iosys('read integer "optimization cycle" from rwf',1,
     $              bscycl,0,' ')
         call iosys('read character search_type from rwf',-1,0,0,srcht)
         call iosys('read character source_d2e from rwf',-1,0,0,srcd2e)
         call iosys('read integer /bsinfo/ from rwf',lenbs,energy,0,' ')
         call iosys('read real energy from rwf',1,energy,0,' ')
         call iosys('read real "bs x" from rwf',nvar,x,0,' ')
         call iosys('read real forces from rwf',nvar,f,0,' ')
         call iosys('read real "bs fc" from rwf',2*nvv,fc,0,' ')
         call iosys('read real force_constants from rwf',nvv,frcnst,
     $              0,' ')
         call iosys('read real "bs xx" from rwf',nvar*mxpt,xx,0,' ')
         call iosys('read real "bs ff" from rwf',nvar*mxpt,ff,0,' ')
         call iosys('read real "bs fs" from rwf',mxpt,fs,0,' ')
         call iosys('read integer "bs ic" from rwf',mxpt,ic,0,' ')
         nwords=nvar*len(vname(1))
         call iosys('read character "variable names" from rwf',nwords,
     $              0,0,vname)
      endif
c
c
      return
      end
