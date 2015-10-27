*deck @(#)dxgrad.f	5.1  11/6/94
      subroutine dxgrad(nvar,f,xnew,dxrms,dxmax,ok)
c***begin prologue     dxgrad.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)dxgrad.f	5.1   11/6/94
c***purpose            
c***description
c
c     in the event that all else has failed, this routine sets
c     up a step in the direction of the gradient.  this is sure
c     to lead to eventual lowering of the energy.
c
c***references
c
c***routines called
c
c***end prologue       dxgrad.f
      implicit none
c     --- input variables -----
      integer nvar
c     --- input arrays (unmodified) ---
      real*8 f(nvar)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 xnew(nvar)
c     --- output variables ---
      logical ok
      real*8 dxrms,dxmax
c     --- scratch arrays ---
c     --- local variables ---
      real*8 scale
      real*8 sdot,arrmax
      parameter (scale=0.5d00)
c
c
      call smul(xnew,f,scale,nvar)
      dxrms=sqrt(sdot(nvar,xnew,1,xnew,1)/float(nvar))
      dxmax=arrmax(xnew,nvar)
      ok=.true.
c
c
      return
      end
