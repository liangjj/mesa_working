*deck @(#)dfbydx.f	5.1  11/6/94
      subroutine dfbydx(ndeg,nfd,a,ap)
c***begin prologue     dfbydx.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)dfbydx.f	5.1   11/6/94
c***purpose            
c***description
c
c     computes the coefficients of the first derivative of a polynomial
c     stored in a.  the form of the polynomial is f=a(1)*x**ndeg+...
c     the first derivative coefficients are returned in ap.
c
c***references
c
c***routines called
c
c***end prologue       dfbydx.f
      implicit none
c     --- input variables -----
      integer ndeg,nfd
c     --- input arrays (unmodified) ---
      real*8 a(ndeg),ap(ndeg)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
c
c
      nfd=ndeg-1
      if(ndeg.ge.1) then
         do 10 i=1,ndeg
            ap(i)=a(i)*float(ndeg-i+1)
   10    continue
      endif
c
c
      return
      end
