*deck @(#)lnrstp.f	5.1  11/6/94
      subroutine lnrstp(iout,prnt,eold,enew,fold,fnew,xout,ok)
c***begin prologue     lnrstp.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)lnrstp.f	5.1   11/6/94
c***purpose            
c***description
c     fit a polynomial through the last two points and extrapolate
c     to the minimum along the line.  this step is referred to as the
c     linear search.  eold,enew,fold,and fnew are the old and new
c     values and derivatives, respectively.
c
c     the polynomial is
c        a(1)*x**4 +a(2)*x**3 +a(3)*x**2 +a(4)*x
c     where x=-0.5 is the old point with energy eold, and x=0.5 is the
c     new point of energy enew.  the zero of energy is taken as the
c     energy of the midpoint of the line segment defined by the two
c     points.  the coefficients are determined by fitting to the given
c     energies and first derivatives, and requiring that the second
c     derivative just reach zero at its minimum.  the latter constraint
c     implies that a(1)=(3/8)*a(2)**2 /a(3).  the output xout is 0.0
c     for the latest point, =1.0 to step all the way back to the old
c     point.
c
c     if there is no such quartic, a cubic is attempted.  if this
c     does not have a minimum between xout=-2.0 and +1.0, a quadratic
c     is tried.
c
c     for either fit, if the new step is worse(higher) than the old
c     one, a minimum is only accepted if it is in between the two
c     points (between -1.0 and 0.0).
c
c***references
c
c***routines called
c
c***end prologue       lnrstp.f
      implicit none
c     --- input variables -----
      integer iout
      logical prnt
      real*8 eold,enew,fold,fnew,xout
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      logical ok
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer ndeg,i,limi,nda,ndda
      integer is(4,4)
      real*8 a(5),da(5),dda(5),scr1(5),scr2(5),xstart(3)
      real*8 xval(2),yval(4),xmat(4,4),scr(4)
      real*8 zero,one,two,three,four,half,cutoff,xmax,cubmax
      real*8 xkcube,fmid,s,p,ssqm4p,root,xlim,fatx,enew1,emid,eold1
      real*8 x,dnew,dmid,dold,ddnew,ddmid,ddold,valmin,dmin,ddmin
c
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00,three=3.0d+00)
      parameter (half=5.0d-01,four=4.0d+00,cutoff=1.0d-10)
      parameter (xmax=1.0d+05,cubmax=1.0d+00)
      data xstart/-5.0d-01,0.0d+00,5.0d-01/
      save xstart
c
 1000 format(5x,'lnrstp: no minimum found in linear search.')
 1010 format(5x,'lnrstp: eold=',e12.4,' enew=',e12.4,
     $      /5x,'        dold=',e12.4,' dnew=',e12.4)
 1020 format(5x,'lnrstp: a=',5(e12.4))
 1030 format(5x,'lnrstp: values at new, mid, old=',3(e12.4),
     $      /5x,'        derivs at new, mid, old=',3(e12.4),
     $      /5x,'                               =',3(e12.4))
 1050 format(5x,'lnrstp: xstart=',3(e12.4))
 1060 format(5x,'lnrstp:  x,   f at min=',2(e12.4),
     $      /5x,'        df, ddf at min=',2(e12.4))
 1070 format(5x,'lnrstp: ssqm4p is negative ... using cubic minimum.')
 1080 format(5x,'lnrstp: a(3) is too small  ... using cubic minimum.')
 1090 format(5x,'lnrstp: x is unacceptable.')
 1100 format(5x,'lnrstp: cubic fit failed.')
c
c     --- evaluate first derivative components old and new along
c         direction of motion.
      if(prnt) write(iout,1010) eold,enew,fold,fnew
      xkcube=eold-enew
      fmid=(fold+fnew)/two
      a(2)=two*(xkcube+fmid)
      a(4)=-(three*xkcube+fmid)/two
      s=fnew-fold
c
c     --- do not attempt extrapolation if a negative second derivative.
c         is encountered.
      ok=.false.
      if(s.ge.cutoff) then
c        --- extrapolate using a quartic curve fitted to first three
c            derivatives and adjusted so that the second derivative
c            just reaches zero but does not become negative.
         p=(three/four)*a(2)**2
         ssqm4p=s*s-four*p
         if(ssqm4p.ge.0) then
            root=sqrt(ssqm4p)
            a(3)=(s+root)/four
            if(a(3).ge.cutoff) then
               a(1)=(s-root)*half
               a(5)=zero
               ndeg=4
               xlim=xmax+half
            else
               if(prnt) write(iout,1080)
            endif
         else
            if(prnt) write(iout,1070)
         endif
c
c        --- if the cubic term is too large, or the step becomes too large,
c            try a straight cubic.
         if(ssqm4p.lt.zero.or.a(3).lt.cutoff) then
            xval(1)=-half
            xval(2)=half
            yval(1)=eold
            yval(2)=enew
            yval(3)=fold
            yval(4)=fnew
            call polfit(2,2,xval,yval,xmat,is,scr,a,ok)
            if(ok) then
               ndeg=3
               xlim=cubmax+half
            else
               if(prnt) write(iout,1100)
               return
            endif
         endif
c
c        --- find the minimum of the polynomial.
         limi=ndeg+1
         if(prnt) write(iout,1020) (a(i),i=1,limi)
         enew1=fatx(ndeg,half,a)
         emid=fatx(ndeg,zero,a)
         eold1=fatx(ndeg,-half,a)
         call dfbydx(ndeg,nda,a,da)
         dnew=fatx(nda,half,da)
         dmid=fatx(nda,zero,da)
         dold=fatx(nda,-half,da)
         call dfbydx(nda,ndda,da,dda)
         ddnew=fatx(ndda,half,dda)
         ddmid=fatx(ndda,zero,dda)
         ddold=fatx(ndda,-half,dda)
         if(prnt) then
            write(iout,1030) enew1,emid,eold1,dnew,dmid,dold,
     $                       ddnew,ddmid,ddold
            write(iout,1050) (xstart(i),i=1,3)
         endif
         call polmin(ndeg,a,3,xstart,da,dda,scr1,scr2,x,ok)
         if(prnt.and.(.not.ok)) write(iout,1000)
         ok=ok.and.
     $     ((abs(x).le.xlim).and.(enew.lt.eold.or.abs(x).le.half))
c
c        --- compute the step vector.
         if(ok) then
            xout=x-half
            valmin=fatx(ndeg,x,a)
            dmin=fatx(nda,x,da)
            ddmin=fatx(ndda,x,dda)
            if(prnt) write(iout,1060) x,valmin,dmin,ddmin
         else
            if(prnt) write(iout,1090)
         endif
      endif
c
c
      return
      end
