*deck @(#)msmain.f	5.2  11/28/95
      subroutine msmain(nvar,nvv,nz,mxpt,toang,ops,x,f,xx,ff,fc,frcnst,
     $                  fs,vname,ic,hinv,lbl,convrg,abnrml,z,zsq,a)
c***begin prologue     msmain.f
c***date written       yymmdd  
c***revision date      11/28/95      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)msmain.f	5.2   11/28/95
c***purpose            
c***description
c
c     local arrays of note.
c     x      ... coordinate vector; real(nvar).
c     f      ... forces(-de/dx); real(nvar).
c     frcnst ... second derivative matrix if evaluated externally.
c                real(nvv).
c     vname  ... variable names vector; character*16(nvar).
c     fc     ... second derivative matrices; real(2*nvv).
c     fs     ... value of the function at earlier points; real(mxpt).
c     ic     ... cycle pointers associated with fs; integer(mxpt).
c     xx     ... previous coordinates; real(nvar,mxpt)
c     ff     ... previous forces; real(nvar,mxpt)
c     hinv   ... inverse hessian; real(nvar,nvar).
c
c     common/msinf/ energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,
c    $              eigmax,eigmin,fncerr,grderr,fnccnv,alpha,
c    $              mxcycl,mscycl,np,neg,
c    $              prnt,tstcrv,updrwf,chkpt,clnup
c     energy ... the value of the function at the current point.
c     convf  ... convergence criterion.
c     fmaxt  ... maximum force allowed before aborting.
c     dxmaxt ... maximum allowed displacement before rescaling.
c     rmax   ... see gshmdt, currently not used.
c     rmin   ...
c     rlim   ...
c     eigmax ... maximum allowed force constant eigenvalue before abort.
c     eigmin ... minimum allowed force constant eigenvalue before abort.
c     fncerr ... estimate of error in value of function.
c     grderr ... estimate of error in value of gradient.
c     fnccnv ... unit conversion factor for function.
c     mxcycl ... maximum number of cycles allowed.
c     mscycl ... current cycle number.
c     np     ... number of points stored.
c     neg    ... order of stationary point for which to search.
c     prnt   ... print switch.
c     tstcrv ... test stationary point for proper curvature?
c     updrwf ... update the rwf?
c     chkpt  ... checkpoint the run?
c     clnup  ... clean the d2e files on the rwf?
c
c***references
c
c***routines called
c
c***end prologue       msmain.f
      implicit none
c     --- input variables -----
      integer nvar,nvv,nz,mxpt
      real*8 toang
      character*(*) ops
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
c     --- input arrays (scratch) ---
      integer ic(mxpt),lbl(nz),a(7*nz+2*nvar)
      real*8 x(nvar),f(nvar),frcnst(nvv),hinv(nvar,nvar)
      real*8 xx(nvar,mxpt),ff(nvar,mxpt),fc(2*nvv),fs(mxpt)
      real*8 z(6*nvar),zsq(nvar*nvar)
c     --- output arrays ---
c     --- output variables ---
      logical convrg,abnrml
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer mxcycl,mscycl,np,neg,ipsav1,ipsav2
      integer intkey
      character*8 srcd2e,namchk*128
      character*4 ians
      logical prnt,exit,tstcrv
      logical updrwf,chkpt,clnup,logkey
      real*8 fpkey
      real*8 convf,dxmaxt,fmaxt,eigmax,eigmin
      real*8 energy,rmax,rmin,rlim,fncerr,grderr,fnccnv,alpha
c
      common/io/inp,iout
      common/msinf/ energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,
     $              eigmax,eigmin,fncerr,grderr,fnccnv,alpha,
     $              mxcycl,mscycl,np,neg,ipsav1,ipsav2,
     $              prnt,tstcrv,updrwf,chkpt,clnup
c
 1000 format(1x,'murtaugh-sargent optimization:')
 1020 format(5x,'restart from: ',a8)
c
c     --- announce our presence
      write(iout,1000)
      exit=.false.
      call iosys('read integer zlbl from rwf',nz,lbl,0,' ')
c
c     --- determine the type of entry.  '"optimization cycle"' has not
c         been written to the rwf is this is a restart or the first point
c         of a new optimization.
      call iosys('does "optimization cycle" exist on rwf',0,0,0,ians)
      if (ians.eq.'no') then
         if(logkey(ops,'opt=restart',.false.,' ')) then
c           --- initial entry, restart optimization.
            call iosys('read character "checkpoint filename" from rwf',
     $           -1,0,0,namchk)
            write(iout,1020) namchk
            call chkpnt('restore')
            call msinfo('read',nvar,nvv,mxpt,vname,x,f,fc,frcnst,xx,
     $                  ff,fs,ic,srcd2e,hinv)
            mxcycl=intkey(ops,'opt=cycles',mxcycl,' ')
            convf=fpkey(ops,'opt=convf',convf,' ')
            dxmaxt=fpkey(ops,'opt=dxmaxt',dxmaxt,' ')
            fmaxt=fpkey(ops,'opt=fmaxt',1.0d+00,' ')
            eigmax=fpkey(ops,'opt=eigmax',eigmax,' ')
            eigmin=fpkey(ops,'opt=eigmin',eigmin,' ')
         else
c           --- initial entry, new optimization.
            call msinit(nvar,nvv,ops,zsq,fc,frcnst,a(7*nz+1),
     $                  x,vname,hinv,nz,lbl,toang,srcd2e,a,a(nz+1),
     $                  a(5*nz+1),a(6*nz+1),z(1),a(7*nz+nvar+1))
            call msinfo('write',nvar,nvv,mxpt,vname,x,f,fc,frcnst,xx,
     $                  ff,fs,ic,srcd2e,hinv)
            convrg=.false.
            return
         endif
      else
c        --- a continuation of the current run.
         call iosys('read integer "optimization cycle" from rwf',
     $              1,mscycl,0,' ')
         call iosys('write character "print flag" to rwf',0,0,0,
     $              'minimum')
         call msinfo('read',nvar,nvv,mxpt,vname,x,f,fc,frcnst,xx,
     $               ff,fs,ic,srcd2e,hinv)
         if(srcd2e.ne.'analytic') then
            call vmove(frcnst,fc,nvv)
            call iosys('write real force_constants to rwf',nvv,frcnst,
     $                 0,' ')
         endif
         if(chkpt) call chkpnt('save')
      endif
c
c     --- if this is the last point at which analytic second derivatives
c         are available, clean up the rwf.
      if(clnup) call clnd2e
c
      call msopt(nz,nvar,nvv,mxpt,toang,vname,x,f,fc,frcnst,xx,
     $            ff,fs,ic,hinv,lbl,convrg,exit,srcd2e,z,z(nvar+1),
     $            zsq,a)
c
c     --- save the current force constant matrix on the chk file.  if the
c         force constant matrix resulted from a frequency run, then it already
c         exists in frcnst.  note that the copy in fc may have been altered by
c         msopt.  if this is a geometry optimization, the matrix fc will have
c         to be copied into frcnst.
      abnrml=exit.and.(.not.convrg)
      if(exit) then
         if(updrwf) call vmove(frcnst,fc,nvv)
      else
c        --- update the rwf and exit the link.
         if(updrwf)  call iosys('write real zvalues to rwf',
     $                           nvar,x,0,' ')
         call msinfo('write',nvar,nvv,mxpt,vname,x,f,fc,frcnst,xx,
     $               ff,fs,ic,srcd2e,hinv)
      endif
c
c
      return
      end
