*deck @(#)bslinr.f	5.1  11/6/94
      subroutine bslinr(nvar,nvv,mxpt,f,xnew,dxrms,dxmax,ok,xout,
     $                  ff,xx,fs,fc,deltax,hdold,hdnew,h)
c***begin prologue     bslinr.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)bslinr.f	5.1   11/6/94
c***purpose            performs linear search for stationary point.
c***description
c     fit a polynomial through the last two points and extrapolate
c     to the minimum along the line.  this step is referred to as the
c     linear search.  the quadratic search is then made starting from
c     this extrapolated minimum in bsquad.  the output xold is 0.0
c     for the latest point with energy fs(1), and -1.0 to step all the
c     way back to the old point with energy fs(2).  neg is the requested
c     curvature.  ok is returned .true. if a step was generated.
c
c     there are three cases in which this routine can perform the
c     linear search.
c
c     1.  neg=0,  ipsav1 or ipsav2=0.  in this case the search is for
c         a local minimum and there are only first derivatives available
c         at one or both of the points.  lnrstp is used to fit a cubic
c         or constrained quartic to locate the minimum.  unless there is
c         only a local maximum between the two points, the quartic will
c         have exactly one minimum.  a minimum in the cubic is only
c         considered valid if it is "reasonably" close to the search
c         region.
c
c
c     2.  neg>0,  ipsav1,ipsav2>0.  here the search is for a higher
c         order stationary point.  the linear search attempts to
c         minimize the norm of the gradient, and consequently second
c         derivatives must be available at both points.  a quartic fit
c         to the gradient norm analogous to method 1 is performed.
c         this is inherently inferior to the newton-raphson step and is
c         only performed if the previous step failed to reduce the norm
c         of the gradient.
c
c     3.  neg=0,  ipsav1,ipsav2>0.  here both first and second
c         derivatives are available at both points.  a general quintic
c         polynomial is fit.  if this does not have a reasonable minimum,
c         method 1 is used instead.
c
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       bslinr.f
      implicit none
c     --- input variables -----
      integer nvar,nvv,mxpt
c     --- input arrays (unmodified) ---
      real*8 ff(nvar,mxpt),xx(nvar,mxpt),fs(mxpt)
c     --- input arrays (scratch) ---
      real*8 deltax(nvar),hdold(nvar),hdnew(nvar),h(nvar,nvar)
      real*8 f(nvar)
      real*8 fc(2*nvv)
c     --- output arrays ---
      real*8 xnew(nvar)
c     --- output variables ---
      logical ok
      real*8 dxrms,dxmax,xout
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer mxcycl,bscycl,np,neg,ipsav1,ipsav2
      integer is(4,6)
      integer i
      logical prnt,tstcrv,updrwf,chkpt,clnup
      real*8 xval(2),yval(6),xmat(6,6),scr(7)
      real*8 a(6),da(6),dda(6),xloc(7),xstart(7)
      real*8 energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,eigmax,eigmin
      real*8 fswtch,fncerr,grderr,fnccnv
      real*8 fold,fnew,gnold,gnnew,fxold,fatx,fxnew,fxmin,dfxold,dfxnew
      real*8 dgnold,dgnnew,dfxmin,ddfxo,ddfxn,ddfxm,rx
      real*8 zero,one
      real*8 sdot,sqrt,arrmax
c
      parameter (zero=0.0d+00,one=1.0d+00)
c
      common/bsinf/ energy,convf,fmaxt,dxmaxt,rmax,rmin,rlim,
     $              eigmax,eigmin,fswtch,fncerr,grderr,fnccnv,
     $              mxcycl,bscycl,np,neg,ipsav1,ipsav2,
     $              prnt,tstcrv,updrwf,chkpt,clnup
      common/io/inp,iout
c
      data xstart/-2.0d0,-1.0d0,-0.75d0,-0.50d0,-0.25d0,0.0d0,1.0d0/
      save xstart
c
 1000 format(5x,'bslinr:  eold= ',e12.4,' enew= ',e12.4,/,
     $       5x,'         fold= ',e12.4,' fnew= ',e12.4,/,
     $       5x,'        gnold= ',e12.4,'gnnew= ',e12.4)
 1010 format(5x,'bslinr:  minimizing gradient -- ipsav =',2i4)
 1020 format(5x,'bslinr:  xval = ',2(e12.4),/,
     $       5x,'         yval = ',3(e12.4),/,
     $       5x,'                ',3(e12.4))
 1030 format(5x,'bslinr:  quintic fit failed.')
 1040 format(5x,'bslinr:  quintic coefficients= ',3(e12.4),
     $      /5x,'                               ',3(e12.4),
     $      /5x,'                       values= ',3(e12.4),
     $      /5x,'                               ',3(e12.4),
     $      /5x,'                               ',3(e12.4))
 1050 format(5x,'bslinr:  no minimum found in quintic.')
 1060 format(5x,'bslinr:  the minimum in the quintic is at ',e12.4)
 1070 format(5x,'bslinr:  using constrained quartic.')
 1081 format(5x,'bslinr:  fold = ',4d12.4,
     $      /5x,'                ',4d12.4)
 1082 format(5x,'bslinr:  fnew = ',4d12.4,
     $      /5x,'                ',4d12.4)
 1083 format(5x,'bslinr:  delx = ',4d12.4,
     $      /5x,'                ',4d12.4)
 1084 format(5x,'bslinr:  hdo  = ',4d12.4,
     $      /5x,'                ',4d12.4)
 1085 format(5x,'bslinr:  hdn  = ',4d12.4,
     $      /5x,'                ',4d12.4)
 1090 format(5x,'bslinr:  hessian for point ',i3)
c
c     --- evaluate first derivative components old and new along
c         direction of motion.
      ok=.false.
      call vsub(deltax,xx(1,2),xx(1,1),nvar)
      fold=sdot(nvar,ff(1,2),1,deltax,1)
      fnew=sdot(nvar,ff(1,1),1,deltax,1)
      gnold=sdot(nvar,ff(1,2),1,ff(1,2),1)
      gnnew=sdot(nvar,ff(1,1),1,ff(1,1),1)
      if(prnt) then
         write(iout,1000) fs(2),fs(1),fold,fnew,gnold,gnnew
         write(iout,1081) (ff(i,2),i=1,nvar)
         write(iout,1082) (ff(i,1),i=1,nvar)
         write(iout,1083) (deltax(i),i=1,nvar)
      endif
      if(ipsav2.ne.0) then
         call trtosq(h,fc(nvv+1),nvar,nvv)
         call ebc(hdold,h,deltax,nvar,nvar,1)
         if(prnt) then
            write(iout,1090) ipsav2
            call matprt(h,nvar,nvar,nvar,nvar,0,0,' ',' ',0,hdold,
     $                  .false.)
            write(iout,1084) (hdold(i),i=1,nvar)
         endif
      endif
      if(ipsav1.ne.0) then
         call trtosq(h,fc(1),nvar,nvv)
         call ebc(hdnew,h,deltax,nvar,nvar,1)
         if(prnt) then
            write(iout,1090) ipsav1
            call matprt(h,nvar,nvar,nvar,nvar,0,0,' ',' ',0,hdnew,
     $                  .false.)
            write(iout,1085) (hdnew(i),i=1,nvar)
         endif
      endif
c
c
      if((ipsav1.eq.0.or.ipsav2.eq.0).and.neg.eq.0) then
c        --- case 1:constrained quartic.
         if(prnt) write(iout,1070)
         call lnrstp(iout,prnt,fs(2)*fnccnv,fs(1)*fnccnv,fold,fnew,
     $               xout,ok)
      else if(neg.ne.0) then
c        --- case 2:minimize the norm of the gradient.
         if(prnt) write(iout,1010) ipsav1,ipsav2
         if(ipsav1.ne.0.and.ipsav2.ne.0) then
            dgnold=sdot(nvar,hdold,1,ff(1,2),1)
            dgnnew=sdot(nvar,hdnew,1,ff(1,1),1)
            call lnrstp(iout,prnt,gnold,gnnew,dgnold,dgnnew,xout,ok)
         endif
      else
c        --- case 3:quintic fit.  if this fails, do a constrained quartic.
         xval(1) = -one
         xval(2) = zero
         yval(1) = fs(2)*fnccnv
         yval(2) = fs(1)*fnccnv
         yval(3) = fold
         yval(4) = fnew
         yval(5) = sdot(nvar,hdold,1,deltax,1)
         yval(6) = sdot(nvar,hdnew,1,deltax,1)
         if(prnt) write(iout,1020) xval,yval
         call polfit(2,3,xval,yval,xmat,is,scr,a,ok)
         if(ok) then
            call polmin(5,a,7,xstart,da,dda,scr,xloc,xout,ok)
            if(.not.ok) xout=zero
            fxold  =fatx(5,-one,a)
            fxnew  =fatx(5,zero,a)
            fxmin  =fatx(5,xout,a)
            dfxold =fatx(4,-one,da)
            dfxnew =fatx(4,zero,da)
            dfxmin =fatx(4,xout,da)
            ddfxo  =fatx(3,-one,dda)
            ddfxn  =fatx(3,zero,dda)
            ddfxm  =fatx(3,xout,dda)
            if(prnt) then
               write(iout,1040) a,fxold,fxnew,fxmin,dfxold,dfxnew,
     $                          dfxmin,ddfxo,ddfxn,ddfxm
               if(ok) then
                  write(iout,1060) xout
               else
                  write(iout,1050)
               endif
            endif
            ok=ok.and.
     $         (fs(1).lt.fs(2).or.(xout.ge.-one.and.xout.le.zero))
         endif
         if(.not.ok) then
c           --- quintic fit didn't work.  try constrained quartic.
            if(prnt) write(iout,1070)
            call lnrstp(iout,prnt,fs(2)*fnccnv,fs(1)*fnccnv,fold,fnew,
     $                  xout,ok)
         endif
      endif
c
c     --- the step xout has been determined. compute xnew.
      if(ok) then
         rx=sdot(nvar,deltax,1,deltax,1)
         call vwxs(xnew,xnew,deltax,-xout,+1,nvar)
         dxrms=sqrt(xout*xout*sdot(nvar,deltax,1,deltax,1)/float(nvar))
         dxmax=abs(xout)*arrmax(deltax,nvar)
         call vwxs(f,f,ff(1,2),-xout,+1,nvar)
         call vwxs(f,f,ff(1,1),xout,+1,nvar)
         fnew=sdot(nvar,f,1,deltax,1)/rx
         call vwxs(f,f,deltax,-fnew,+1,nvar)
      endif
c
c
      return
      end
