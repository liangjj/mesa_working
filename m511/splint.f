*deck @(#)splint.f	5.1 11/6/94
      subroutine splint(xa,ya,y2a,n,x,y,ind,ng)
c***begin prologue     splint.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             
c***source             @(#)splint.f	5.1   11/6/94
c***purpose            
c***description
c   given the arrays xa and ya of length n -- which tabulate a
c   function with the xa(i) in order --
c   and given the array y2a, which is the output from spline3,
c   and given a value of x, find the cubic spline interpolated 
c   value y.  
c   ind points to the nearest knot in the sequence to the left of x.
c
c***references
c
c***routines called
c
c***end prologue       splint.f
      implicit none
c     --- input variables -----
      integer n,ng
      real*8 x(ng)
c     --- input arrays (unmodified) ---
      integer ind(ng)
      real*8 xa(n),ya(n),y2a(n)      
c     --- input arrays (scratch) ---
c     --- output arrays ---
c     --- output variables ---
      real*8 y(ng)
c     --- scratch arrays ---
c     --- local variables ---
      real*8 a,b,h,zero,one,six
      integer i
c
      parameter (zero=0.0d+00,one=1.0d+00,six=6.0d+00)
c
c     --- get the step length.
      do 10 i=1,ng
         h=xa(ind(i)+1)-xa(ind(i))
         if(h.eq.zero) then
            call lnkerr('spline points must be distinct')
         endif
c
c        --- evaluate the cubic spline.
         a=(xa(ind(i)+1)-x(i))/h
         b=(x(i)-xa(ind(i)))/h
         y(i)=a*ya(ind(i))+b*ya(ind(i)+1)+
     $       ((a**3-a)*y2a(ind(i))+(b**3-b)*y2a(ind(i)+1))*(h*h)/six
   10 continue
c
c
      return
      end
