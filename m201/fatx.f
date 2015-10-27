*deck @(#)fatx.f	5.1  11/6/94
      function fatx(n,x,a)
c***begin prologue     fatx.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)fatx.f	5.1   11/6/94
c***purpose            
c***description
c   this function evaluates the polynomial of degree n with
c   coefficients a at the point x.  the form of the polynomial is
c   f=a(1)*x**n +a(2)*x**(n-1) +... +a(n)*x +a(n+1).
c
c   note that a is used to n+1 locations.
c
c***references
c
c***routines called
c
c***end prologue       fatx.f
      implicit none
      real*8 fatx
c     --- input variables -----
      integer n
c     --- input arrays (unmodified) ---
      real*8 a(n+1)
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
      real*8 x,xx,zero,sum
      parameter (zero=0.0d+00)
c
c     --- allow for a zero order polynomial.
      if(n.lt.0) then
         fatx=zero
      else if(n.eq.0) then
         fatx=a(1)
      else
         sum=a(n+1)
         xx=x
         do 40 i=1,n
            sum=sum+a(n-i+1)*xx
            xx=xx*x
   40    continue
         fatx=sum
      endif
c
c
      return
      end
