*deck @(#)polmin.f	5.1  11/6/94
      subroutine polmin(n,a,npt,xstart,da,dda,val,xloc,xmin,ok)
c***begin prologue     polmin.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)polmin.f	5.1   11/6/94
c***purpose            
c***description
c     this subroutine looks for minima in the polynomial of degree n
c     given by
c        f=a(1)*x**n +.... +a(n+1).
c
c     npt starting points are provided in array xstart. da and dda are
c     loaded with the first and second derivatives of a and should be
c     dimensioned n and n-1, respectively.  val and xloc are scratch
c     arrays dimensioned npt.  the location of the minimum is
c     returned in xmin.
c     ok is returned .true. if a minimum was found.
c
c***references
c
c***routines called
c
c***end prologue       polmin.f
      implicit none
c     --- input variables -----
      integer n,npt
c     --- input arrays (unmodified) ---
c     --- input arrays (scratch) ---
      real*8 a(n+1),xstart(npt),da(n),dda(n-1)
      real*8 val(npt),xloc(npt)
c     --- output arrays ---
c     --- output variables ---
      logical ok
      real*8 xmin
c     --- scratch arrays ---
c     --- local variables ---
      integer nmin,i,istep,nd,ndd,ipt,mxstep
      real*8 x,dx,ddx,cutoff,fval,fatx
      parameter (cutoff=1.0d-10,mxstep=1000)
c
c     --- differentiate the polynomial twice.
      call dfbydx(n,nd,a,da)
      call dfbydx(nd,ndd,da,dda)
c
c     --- for each starting point at which the curvature is reasonable, use
c         newton-raphson to locate the nearest stationary point.
      nmin=0
      ok=.false.
      do 30 ipt=1,npt
         x=xstart(ipt)
         do 10 istep=1,mxstep
            dx=fatx(nd,x,da)
            ddx=fatx(ndd,x,dda)
            if(abs(dx).lt.cutoff) goto 20
            if(abs(ddx).lt.cutoff) goto 30
            x=x-dx/ddx
   10    continue
         return
c
c
   20    if(ddx.ge.cutoff) then
            nmin=nmin+1
            xloc(nmin)=x
            val(nmin)=fatx(n,x,a)
         endif
   30 continue
c
c
      if(nmin.gt.0) then
         ok=.true.
         xmin=xloc(1)
         if(nmin.gt.1) then
            fval=val(1)
            do 40 i=2,nmin
               if(val(i).gt.fval) then
                  fval=val(i)
                  xmin=xloc(i)
               endif
   40       continue
         endif
      endif
c
c
      return
      end
