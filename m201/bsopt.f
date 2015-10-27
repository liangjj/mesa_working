*deck @(#)bsopt.f	5.1  11/6/94
      subroutine bsopt(nz,nvar,nvv,mxpt,toang,vname,x,f,fc,frcnst,xx,
     $                  ff,fs,ic,lbl,convrg,exit,srcd2e,srcht,xnew,
     $                  xquad,ftemp,z,zsq,a)
c***begin prologue     bsopt.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al.(g82)
c***source             @(#)bsopt.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       bsopt.f
      implicit none
c     --- input variables -----
      integer nz,nvar,nvv,mxpt
      character*(*) srcd2e,srcht
      real*8 toang
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 x(nvar),f(nvar)
c     --- output variables ---
      logical convrg,exit
c     --- scratch arrays ---
      integer ic(mxpt),lbl(nz),a(nvar)
      real*8 fc(2*nvv),frcnst(nvv),xx(nvar,mxpt)
      real*8 ff(nvar,mxpt),fs(mxpt)
      real*8 xnew(nvar),xquad(nvar),ftemp(nvar)
      real*8 z(6*nvar),zsq(2*nvar*nvar)
c     --- local variables ---
      integer inp,iout
      integer nnege,i,j,ij
      integer mxcycl,bscycl,np,neg,ipsav1,ipsav2
      character tname*16
      character result*3
      logical linmin,lints,linr,didlin,rises,quad,steppd,ok
      logical prnt,tstcrv
      logical updrwf,chkpt,clnup
      real*8 convf,dxmaxt,fswtch,fmaxt,eigmax,eigmin
      real*8 energy,rmax,rmin,rlim,fncerr,grderr,fnccnv
      real*8 zero,two,cnvfmx,cnvfx
      real*8 sdot,sqrt,arrmax,fmax,frms,dxmax,dxrms
      real*8 xlin,scale,de,convx,t1,t2
c
      common/bsinf/ energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,
     $              eigmax,eigmin,fswtch,fncerr,grderr,fnccnv,
     $              mxcycl,bscycl,np,neg,ipsav1,ipsav2,
     $              prnt,tstcrv,updrwf,chkpt,clnup
      common/io/inp,iout
c
      parameter (zero=0.0d+00,two=2.0d+00,cnvfmx=1.5d+00,cnvfx=4.0d+00)
c
 1000 format(5x,'convergence tests:',
     $      /10x,16x,'item',8x,'value',4x,'threshold',2x,'converged?')
 1002 format(10x,7x,'maximum force',5x,f9.6,5x,f9.6,5x,a3)
 1004 format(10x,7x,'    rms force',5x,f9.6,5x,f9.6,5x,a3)
 1006 format(10x,'maximum displacement',5x,f9.6,5x,f9.6,5x,a3)
 1008 format(10x,'    rms displacement',5x,f9.6,5x,f9.6,5x,a3)
 1010 format(5x,'the second derivative matrix:')
 1015 format(5x,'energy',12x,f20.6,/5x,'predicted change',2x,f20.6)
 1016 format(5x,'search for a local minimum.')
 1017 format(5x,'search for a saddle point.')
 1018 format(5x,'search for a stationary point of order ',i2,'.')
 1020 format(5x,'cycle number ',i3,' out of a maximum of ',i3,
     $      /5x,'all quantities printed in hartrees-bohr-radians.')
 1030 format(5x,'step:')
 1050 format(11x,'variable',7x,'old x',4x,'-de/dx',3(3x,'delta x'),5x,
     $        'new x',
     $      /41x,'(linear)',4x,'(quad)',3x,'(total)',
     $      /)
 1060 format(3x,a16,2x,6f10.5)
 1080 format(5x,'optimization completed.',
     $      /5x,'-- stationary point found.')
 1090 format(5x,'optimization stopped.'
     $      /5x,'-- number of cycles exceeded,  mxcycl=',i4,
     $      /5x,'flag reset to prevent archiving.')
 1100 format(5x,'optimization stopped.'
     $      /5x,'-- wrong number of negative eigenvalues:',
     $      /5x,' desired=',i3,5x,'actual=',i3,
     $      /5x,'flag reset to prevent archiving.')
 1130 format(5x,'maximum step size (dxmaxt=',f8.3,
     $          ') exceeded in linear search.'
     $      /5x,'-- step size scaled by ',f7.3,
     $      /5x,'-- skip quadratic or steepest descent search.')
 1200 format(5x,'second derivative matrix not updated -- first cycle.')
 1210 format(5x,'second derivative matrix not updated -- analytic',
     $          ' derivatives used.')
 1220 format(5x,'linear search not attempted -- srcht set to quad.')
 1222 format(5x,'linear search not attempted -- no second derivatives',
     $          ' available.')
 1224 format(5x,'linear search not attempted -- first point.')
 1226 format(5x,'linear search not attempted -- rms force is less than',
     $          ' fswtch (',f8.5,').')
 1227 format(5x,'linear search not attempted -- energy rises or',
     $          ' forces went down.')
 1230 format(5x,'energy rises -- skip quadratic search.')
 1235 format(5x,'energy rises -- skip steepest descent search.')
 1237 format(5x,'gradient rises -- skip quadratic search.')
 1240 format(5x,'maximum step size (dxmaxt=',f8.3,
     $          ') exceeded in quadratic search.'
     $      /5x,'-- step size scaled by ',f7.3)
 2010 format(5x,'optimization aborted.',
     $      /5x,'-- gradient out of range.',
     $      /5x,'-- maximum allowed force (fmaxt=',f8.3,').')
 2020 format(5x,'optimization aborted.',
     $      /5x,'-- no acceptable step.')
 2030 format(5x,'linear search skipped for unknown reason.')
 2040 format(5x,'the linear search produced a scale factor of',f8.5,'.')
 2050 format(5x,'skip linear search -- no minimum in search direction.')
 2070 format(5x,'steepest descent instead of quadratic search.')
c
c     --- increment the cycle counter and print topology information.
      bscycl=bscycl+1
      if(neg.eq.0) then
         write(iout,1016)
      else if(neg.eq.1) then
         write(iout,1017)
      else if(neg.gt.1) then
         write(iout,1018) neg
      endif
      write(iout,1020) bscycl,mxcycl
c
c     --- set the number of negative eigenvalues to be accepted in the search,
c         compute the rms force, and test convergence.
      nnege=neg
      fmax=arrmax(f,nvar)
      frms=sqrt(sdot(nvar,f,1,f,1)/float(nvar))
      exit=(frms.lt.convf).and.(fmax.lt.cnvfmx*convf)
      if(fmax.ge.fmaxt.and.updrwf) then
         call matprt(f,1,nvar,1,nvar,1,1,'-de/dx',vname,0,xnew,
     $               .false.)
         write(iout,2010) fmaxt
         call lnkerr(' ')
      endif
c
c     --- prepare arrays for search, and save the present point.
      call rzero(xnew,nvar)
      call rzero(xquad,nvar)
      call vmove(ftemp,f,nvar)
      call savept(nvar,nvv,mxpt,rises,ff,xx,f,x,fc,frcnst,ic,
     $            fs,z,z(nvar+1),srcd2e,np,bscycl,neg,
     $            ipsav1,ipsav2,energy,fncerr)
c
c     --- update the second derivatives.
      if(srcd2e.ne.'analytic') then
         if(np.eq.1) then
            write(iout,1200)
         else if(np.gt.1) then
            call updtd2(nvar,nvv,mxpt,xx,x,ff,f,fc,fs,ic,z,
     $                  z(nvar+1),zsq,a,ipsav1,np,rmax,rmin,rlim,grderr)
         endif
      else
         write(iout,1210)
      endif
c
c     --- print the second derivative matrix.
      write(iout,1010)
      call trtosq(zsq,fc,nvar,nvv)
      call matprt(zsq,nvar,nvar,nvar,nvar,1,1,vname,vname,1,fc,
     $            .false.)
c
c     --- determine the next step.
c         a linear step is taken only if we are allowed to do it,
c         at least two points have been computed, and either the search is for
c         a local minimum or there are analytic second derivatives at both
c         points and the gradient has risen.
      steppd=.false.
      linmin=(neg.eq.0).and.(srcht.eq.'lq'.or.srcht.eq.'lsd'
     $                       .or.(srcht.eq.'lsd?'.and.rises))
      lints =(neg.ne.0).and.(srcht.eq.'lq'.and.rises
     $                       .and.ipsav1.ne.0.and.ipsav2.ne.0)
      linr  =(np.gt.1).and.(frms.gt.fswtch).and.(linmin.or.lints)
      quad  =index(srcht,'q').ne.0
c
c     --- get the step.
      if(linr) then
         call bslinr(nvar,nvv,mxpt,ftemp,xnew,dxrms,dxmax,didlin,xlin,
     $               ff,xx,fs,fc,z,z(nvar+1),z(2*nvar+1),
     $               zsq)
         steppd=steppd.or.didlin
         if(dxmax.gt.dxmaxt) then
            scale=dxmaxt/dxmax
            xlin=xlin*scale
            write(iout,1130) dxmaxt,scale
            call smul(xnew,xnew,scale,nvar)
            dxmax=dxmaxt
         else
            dxmax=zero
            if(rises) then
               if(lints) write(iout,1237)
               if(quad) then
                  if(linmin) write(iout,1230)
               else
                  if(neg.eq.0) write(iout,1235)
               endif
            endif
            if(quad.and.(.not.rises.or.(.not.didlin.and.neg.ne.0))) then
               call bsquad(nvar,nvv,ftemp,xquad,dxrms,dxmax,ok,nnege,
     $                     fc,vname,z,zsq,zsq(nvar*nvar+1),z(nvar+1))
               steppd=steppd.or.ok
            endif
         endif
      else
         call bsquad(nvar,nvv,ftemp,xquad,dxrms,dxmax,ok,nnege,
     $               fc,vname,z,zsq,zsq(nvar*nvar+1),z(nvar+1))
         steppd=steppd.or.ok
c
c        --- print the reason for skipping the linear search here as bsquad
c            prints the eigenvalues of the second derivative matrix.
         if(srcht.eq.'q') then
            write(iout,1220)
         else if(np.le.1) then
            write(iout,1224)
         else if(rises.and.neg.ne.0.and.(ipsav1.eq.0.or.ipsav2.eq.0))
     $           then
            write(iout,1222)
         else if(.not.rises.and.(srcht.eq.'lsd'.or.neg.ne.0))then
            write(iout,1227)
         else if(frms.lt.fswtch) then
            write(iout,1226) fswtch
         else
            write(iout,2030)
            call lnkerr(' ')
         endif
      endif
c
c     --- if both the linear and quadratic steps have failed, or we're doing
c         a steepest descent search, step in the direction of the gradient.
      if(linr) then
         if(didlin) then
            write(iout,2040) xlin
         else
            write(iout,2050)
         endif
      endif
      if(neg.eq.0.and.(.not.steppd.or.(.not.rises.and..not.quad))) then
         write(iout,2070)
         call dxgrad(nvar,f,xquad,dxrms,dxmax,ok)
         steppd=steppd.or.ok
      endif
      if(dxmax.gt.dxmaxt) then
         scale=dxmaxt/dxmax
         write(iout,1240) dxmaxt,scale
         call smul(xquad,xquad,scale,nvar)
      endif
c
c     --- test the size of the displacement, update the position vector,
c         and save original values for printing.
      call vmove(z,x,nvar)
      call vmove(z(nvar+1),xnew,nvar)
      call vadd(xnew,xnew,xquad,nvar)
      call vadd(x,x,xnew,nvar)
      dxmax=arrmax(xnew,nvar)
      dxrms=sqrt(sdot(nvar,xnew,1,xnew,1)/float(nvar))
c
c     --- predict the change in energy.
      de=zero
      ij=0
      do 30 i=1,nvar
         do 20 j=1,i
            ij=ij+1
            de=de-xnew(i)*fc(ij)*xnew(j)
   20    continue
         de=de+xnew(i)*fc(ij)*xnew(i)/two
   30 continue
      de=de/fnccnv
c
c     --- test convergence.
      convx=cnvfx*convf
      exit=exit.and.(dxrms.lt.convx).and.(dxmax.lt.cnvfmx*convx)
     $         .and.(nnege.eq.neg.or..not.tstcrv)
c
c     --- output.
      write(iout,1030)
      write(iout,1050)
      do 40 i=1,nvar
         call crjust(vname(i),tname)
         write(iout,1060) tname,z(i),f(i),z(nvar+i),xquad(i),
     $                    xnew(i),x(i)
   40 continue
      t1=convf*cnvfmx
      t2=convx*cnvfmx
      write(iout,1000)
      call convgd(fmax,t1,result)
      write(iout,1002) fmax,t1,result
      call convgd(frms,convf,result)
      write(iout,1004) frms,convf,result
      call convgd(dxmax,t2,result)
      write(iout,1006) dxmax,t2,result
      call convgd(dxrms,convx,result)
      write(iout,1008) dxrms,convx,result
      write(iout,1015) energy,de
c
c     --- if we haven't converged, check for exceeding the maximum number of
c         cycles, the wrong number of negative eigenvalues, or a failure to
c         take a step.
c
      if(exit) then
c        --- if we've converged, print final parameters.
         convrg=.true.
         write(iout,1080)
         do 50 i=1,nvar
            x(i)=x(i)-xnew(i)
            a(i)=99
   50    continue
         call prmtbl(1,vname,x,a,f,nvar,lbl,nz,toang)
      else
c        --- no convergence, shall we proceed?
         convrg=.false.
         if(.not.steppd) then
            write(iout,2020)
            exit=.true.
         else if(nnege.ne.neg.and.tstcrv.and.(.not.updrwf)) then
            write(iout,1100) neg,nnege
            exit=.true.
         else if(bscycl.ge.mxcycl) then
            write(iout,1090) mxcycl
            exit=.true.
         endif
      endif
c
c
      if(exit.and.(.not.convrg)) then
c        --- abnormal exit
         call iosys('write integer archive to rwf',1,.false.,0,' ')
         do 80 i=1,nvar
            a(i)=99
   80    continue
         call prmtbl(2,vname,x,a,f,nvar,lbl,nz,toang)
      else if(.not.convrg) then
c        --- no exit, but print non-optimized parameters each iteration.
         do 90 i=1,nvar
            a(i)=99
   90    continue
         call prmtbl(2,vname,x,a,f,nvar,lbl,nz,toang)
      endif
c
c
      return
      end
