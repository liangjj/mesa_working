*deck @(#)msopt.f	5.1  11/6/94
      subroutine msopt(nz,nvar,nvv,mxpt,toang,vname,x,f,fc,frcnst,xx,
     $                  ff,fs,ic,hinv,lbl,convrg,exit,srcd2e,xnew,
     $                  z,zsq,a)
c***begin prologue     msopt.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)msopt.f	5.1   11/6/94
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       msopt.f
      implicit none
c     --- input variables -----
      integer nz,nvar,nvv,mxpt
      character*(*) srcd2e
      real*8 toang
c     --- input arrays (unmodified) ---
      character*(*) vname(nvar)
c     --- input arrays (scratch) ---
      integer ic(mxpt),lbl(nz),a(nvar)
      real*8 fc(2*nvv),frcnst(nvv),xx(nvar,mxpt)
      real*8 ff(nvar,mxpt),fs(mxpt),hinv(nvar,nvar)
      real*8 xnew(nvar)
      real*8 z(5*nvar),zsq(nvar*nvar)
c     --- output arrays ---
      real*8 x(nvar),f(nvar)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer mxcycl,mscycl,np,neg,ipsav1,ipsav2
      integer nnege,i,j,ij
      character tname*16
      character result*3
      logical rises,updd2e,steppd,goodstp,convrg
      logical prnt,exit,tstcrv
      logical updrwf,chkpt,clnup
      real*8 convf,dxmaxt,fmaxt,eigmax,eigmin
      real*8 energy,rmax,rmin,rlim,fncerr,grderr,fnccnv,alpha,two
      real*8 zero,cnvfmx,cnvfx,convx,one,detmin
      real*8 fmax,frms,det,dxmax,dxrms,scale,de,t1,t2
      real*8 arrmax,sdot,sqrt
c
      parameter (zero=0.0d+00,two=2.0d+00,cnvfmx=1.5d+00,cnvfx=4.0d+00)
      parameter (one=1.0d+00,detmin=1.0d-06)
c
      common/msinf/ energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,
     $              eigmax,eigmin,fncerr,grderr,fnccnv,alpha,
     $              mxcycl,mscycl,np,neg,ipsav1,ipsav2,
     $              prnt,tstcrv,updrwf,chkpt,clnup
      common/io/inp,iout
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
 1019 format(5x,'gradient only -- no optimization.')
 1020 format(5x,'cycle number ',i3,' out of a maximum of ',i3,
     $      /5x,'all quantities printed in hartrees-bohr-radians.')
 1030 format(5x,'step:')
 1050 format(11x,'variable',7x,'old x',4x,'-de/dx',3x,'delta x',5x,
     $        'new x',/)
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
 1200 format(5x,'second derivative matrix not updated -- first cycle.')
 1205 format(5x,'update second derivatives using information from',
     $       ' cycles',5i3,(/5x,53x,5i3))
 1210 format(5x,'analytic second derivatives:take newton-raphson step')
 1215 format(5x,'inverse hessian is singular. det=',e12.4)
 1216 format(5x,'analytic hessian is singular. det=',e12.4)
 1240 format(5x,'maximum step size (dxmaxt=',f8.3,
     $          ') exceeded in search.'
     $      /5x,'-- step size scaled by ',f7.3)
 2010 format(5x,'optimization aborted.',
     $      /5x,'-- gradient out of range.',
     $      /5x,'-- maximum allowed force (fmaxt=',f8.3,').')
 2020 format(5x,'optimization aborted.',
     $      /5x,'-- no acceptable step.',
     $      /5x,'--alpha(',f9.6,') is too small.')
c
c     --- increment the cycle counter and print topology information.
      mscycl=mscycl+1
      if(mxcycl.eq.1) write(iout,1019)
      if(neg.eq.0) then
         if(mxcycl.gt.1) write(iout,1016)
      else if(neg.eq.1) then
         write(iout,1017)
      else if(neg.gt.1) then
         write(iout,1018) neg
      endif
      if(mxcycl.gt.1) write(iout,1020) mscycl,mxcycl
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
c     --- if analytic second derivatives are available, move them into fc,
c         and invert the hessian.
      if(srcd2e.eq.'analytic') then
         call vmove(fc,frcnst,nvv)
         call trtosq(hinv,fc,nvar,nvv)
         call minvrt(hinv,nvar,nvar,det,a,z)
         if(abs(det).le.detmin) write(iout,1216) det
      endif
c
c     --- save the present point, and push the stack.
      if(mscycl.le.1) then
         call rzero(xx(1,1),nvar*mxpt)
         call rzero(ff(1,1),nvar*mxpt)
         call rzero(fs,mxpt)
         call izero(ic,mxpt)
      endif
      call savept(nvar,nvv,mxpt,rises,ff,xx,f,x,fc,frcnst,ic,
     $            fs,z,z(nvar+1),srcd2e,np,mscycl,neg,
     $            ipsav1,ipsav2,energy,fncerr)
c
c     --- prepare for the search.
      call smul(z,ff(1,1),-one,nvar)
      call smul(z(nvar+1),ff(1,2),-one,nvar)
c
c     --- determine the next step.
      call msstep(nvar,mscycl,xnew,xx(1,1),xx(1,2),z,z(nvar+1),hinv,
     $            alpha,z(2*nvar+1),z(3*nvar+1),z(4*nvar+1),fs(1),fs(2),
     $            goodstp,steppd,updd2e,srcd2e,fncerr,iout)
c
c     --- if the current step failed, msstep returns the new step relative to
c         the previous point.  modify the current variables and forces
c         to reflect this.
      if(mscycl.gt.1.and.(.not.goodstp)) then
         call vmove(xx(1,1),xx(1,2),nvar)
         call vmove(ff(1,1),ff(1,2),nvar)
         ic(1)=ic(2)
         call vmove(x,xx(1,2),nvar)
         call vmove(f,ff(1,2),nvar)
      endif
c
c     --- invert the inverse hessian to get updated second derivatives,
c         and print.
      call vmove(zsq,hinv,nvar*nvar)
      call minvrt(zsq,nvar,nvar,det,a,z)
      if(abs(det).le.detmin) then
         write(iout,1215) det
      endif
      call sqtotr(fc,zsq,nvar,nvv)
      if(srcd2e.ne.'analytic') then
         if(np.eq.1) then
            if(mxcycl.gt.1) write(iout,1200)
         else if(np.gt.1) then
            if(updd2e) write(iout,1205) ic(2),ic(1)
         endif
      else
         write(iout,1210)
      endif
      if(mscycl.le.1.and.mxcycl.gt.1) then
         write(iout,1010)
         call matprt(zsq,nvar,nvar,nvar,nvar,1,1,vname,vname,1,fc,
     $               .false.)
      endif
c
c     --- test the size of the displacement, update the position vector,
c         and save original values for printing.
      dxmax=arrmax(xnew,nvar)
      dxrms=sqrt(sdot(nvar,xnew,1,xnew,1)/float(nvar))
      if(dxmax.gt.dxmaxt) then
         scale=dxmaxt/dxmax
         write(iout,1240) dxmaxt,scale
         alpha=alpha*scale
         call smul(xnew,xnew,scale,nvar)
         dxmax=arrmax(xnew,nvar)
         dxrms=sqrt(sdot(nvar,xnew,1,xnew,1)/float(nvar))
      endif
      call vmove(z,x,nvar)
      call vmove(z(nvar+1),xnew,nvar)
      call vadd(x,x,xnew,nvar)
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
c         if we're just calculating the gradient without optimizing,
c         don't bother with this stuff.
      if(mxcycl.gt.1) then
         write(iout,1030)
         write(iout,1050)
      do 40 i=1,nvar
         call crjust(vname(i),tname)
         write(iout,1060) tname,z(i),f(i),z(nvar+i),x(i)
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
      endif
c
c     --- if we haven't converged, check for exceeding the maximum number of
c         cycles, the wrong number of negative eigenvalues, or a failure to
c         take a step.
c         if we've converged, print final parameters.
      if(exit) then
         convrg=.true.
         write(iout,1080)
         do 50 i=1,nvar
            x(i)=x(i)-xnew(i)
            a(i)=99
   50    continue
         call prmtbl(1,vname,x,a,f,nvar,lbl,nz,toang)
      else
         convrg=.false.
         if(.not.steppd) then
            write(iout,2020) alpha
            exit=.true.
         else if(nnege.ne.neg.and.tstcrv.and.(.not.updrwf)) then
            write(iout,1100) neg,nnege
            exit=.true.
         else if(mscycl.ge.mxcycl) then
            if(mxcycl.gt.1) write(iout,1090) mxcycl
            exit=.true.
         endif
      endif
c
c     --- if we are exiting abnormally, print the parameters.
      if(exit.and.(.not.convrg)) then
         call iosys('write integer archive to rwf',1,.false.,0,' ')
         do 80 i=1,nvar
            a(i)=99
   80    continue
         if(mxcycl.gt.1)call prmtbl(2,vname,x,a,f,nvar,lbl,nz,toang)
      else if(.not.convrg) then
c        --- not converged, but print non-optimum parameters each iteration.
         do 90 i=1,nvar
            a(i)=99
   90    continue
         if(mxcycl.gt.1)call prmtbl(2,vname,x,a,f,nvar,lbl,nz,toang)
      endif
c
c
      return
      end
