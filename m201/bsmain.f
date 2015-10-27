*deck @(#)bsmain.f	5.2  4/17/95
      subroutine bsmain(nvar,nvv,nz,mxpt,toang,ops,x,f,xx,ff,fc,frcnst,
     $                  fs,vname,ic,lbl,convrg,abnrml,z,zsq,a)
c***begin prologue     bsmain.f
c***date written       yymmdd  
c***revision date      4/17/95      
c
c***keywords           
c***author             binkley, et al. (g82)j
c***source             @(#)bsmain.f	5.2   4/17/95
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
c
c     common/bsinf/ energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,
c    $              eigmax,eigmin,fswtch,fncerr,grderr,fnccnv,
c    $              mxcycl,bscycl,np,neg,ipsav1,ipsav2,
c    $              prnt,tstcrv,updrwf,chkpt,clnup
c     energy ... the value of the function at the current point.
c     convf  ... convergence criterion.
c     fmaxt  ... maximum force allowed before aborting.
c     dxmaxt ... maximum allowed displacement before rescaling.
c     rmax   ... see gshmdt
c     rmin   ...
c     rlim   ...
c     eigmax ... maximum allowed force constant eigenvalue before abort.
c     eigmin ... minimum allowed force constant eigenvalue before abort.
c     fswtch ... minimum rms force for linear search.
c     fncerr ... estimate of error in value of function.
c     grderr ... estimate of error in value of gradient.
c     fnccnv ... unit conversion factor for function.
c     mxcycl ... maximum number of cycles allowed.
c     bscycl ... current cycle number.
c     np     ... number of points stored.
c     neg    ... order of stationary point for which to search.
c     ipsav1 ... indicates the point for which analytic second
c                derivatives are in the first part of fc.
c                if there are none, it is zero.
c     ipsav2 ... indicates the point for which analytic second
c                derivatives are in the second part of fc.
c                if there are none, it is zero.
c     prnt   ... print switch.
c     tstcrv ... test stationary point for proper curvature?
c     updrwf ... update the rwf?
c     chkpt  ... checkpoint the run?
c     clnup  ... clean the d2e files on the rwf?
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       bsmain.f
      implicit none
c     --- input variables -----
      integer nvar,nvv,nz,mxpt
      real*8 toang
c     --- input arrays (unmodified) ---
      character*(*) ops
c     --- input arrays (scratch) ---
      integer ic(mxpt),lbl(nz),a(7*nz+2*nvar)
      character*(*) vname(nvar)
      real*8 x(nvar),f(nvar),frcnst(nvv)
      real*8 xx(nvar,mxpt),ff(nvar,mxpt),fc(2*nvv),fs(mxpt)
      real*8 z(9*nvar),zsq(2*nvar*nvar)
c     --- output arrays ---
c     --- output variables ---
      logical convrg,abnrml
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer mxcycl,bscycl,np,neg,ipsav1,ipsav2,intkey
      character*128 namchk
      character*8 srcd2e,srcht
      character*4 ians
      logical prnt,exit,tstcrv
      logical updrwf,chkpt,clnup,logkey
      real*8 convf,dxmaxt,fswtch,fmaxt,eigmax,eigmin
      real*8 energy,rmax,rmin,rlim,fncerr,grderr,fnccnv
      real*8 fpkey
c
      common/io/inp,iout
      common/bsinf/ energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,
     $              eigmax,eigmin,fswtch,fncerr,grderr,fnccnv,
     $              mxcycl,bscycl,np,neg,ipsav1,ipsav2,
     $              prnt,tstcrv,updrwf,chkpt,clnup
c
 1000 format(1x,'berny optimization:')
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
            call bsinfo('read',nvar,nvv,mxpt,vname,x,f,fc,frcnst,xx,
     $                  ff,fs,ic,srcd2e,srcht)
            mxcycl=intkey(ops,'opt=cycles',mxcycl,' ')
            fswtch=fpkey(ops,'opt=fswitch',fswtch,' ')
            fmaxt=fpkey(ops,'opt=fmaxt',fmaxt,' ')
         else
c           --- initial entry, new optimization.
            call bsinit(nvar,nvv,ops,zsq,fc,frcnst,a(7*nz+1),
     $                  x,vname,nz,lbl,toang,srcd2e,srcht,a,a(nz+1),
     $                  a(5*nz+1),a(6*nz+1),z(1),a(7*nz+nvar+1))
            call bsinfo('write',nvar,nvv,mxpt,vname,x,f,fc,frcnst,xx,
     $                  ff,fs,ic,srcd2e,srcht)
            convrg=.false.
            return
         endif
      else
c        --- a continuation of the current run.
         call iosys('read integer "optimization cycle" from rwf',
     $        1,bscycl,0,' ')
         call iosys('write character "print flag" to rwf',0,0,0,
     $              'minimum')
         call bsinfo('read',nvar,nvv,mxpt,vname,x,f,fc,frcnst,xx,
     $               ff,fs,ic,srcd2e,srcht)
         if(srcd2e.ne.'analytic') then
            call vmove(frcnst,fc,nvv)
            call iosys('write real force_constants to rwf',nvv,frcnst,
     $              0,' ')
         endif
         if(chkpt) call chkpnt('save')
      endif
c
c     --- if this is the last point at which analytic second derivatives
c         are available, clean up the rwf.
      if(clnup) call clnd2e
c
c     --- call the main driver.
      call bsopt(nz,nvar,nvv,mxpt,toang,vname,x,f,fc,frcnst,xx,
     $           ff,fs,ic,lbl,convrg,exit,srcd2e,srcht,z,z(nvar+1),
     $           z(2*nvar+1),z(3*nvar+1),zsq,a)
c
c     --- save the current force constant matrix on the chk file.  if the
c         force constant matrix resulted from a frequency run, then it already
c         exists in frcnst.  note that the copy in fc may have been altered by
c         bsquad.  if this is a geometry optimization, the matrix fc will have
c         to be copied into frcnst.
      abnrml=exit.and.(.not.convrg)
      if(exit) then
         if(updrwf) call vmove(frcnst,fc,nvv)
      else
c        --- update the rwf and exit the link.
         if(updrwf)  call iosys('write real zvalues to rwf',
     $                           nvar,x,0,' ')
         call bsinfo('write',nvar,nvv,mxpt,vname,x,f,fc,frcnst,xx,
     $               ff,fs,ic,srcd2e,srcht)
      endif
c
c
      return
      end
