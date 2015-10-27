*deck @(#)spline3.f	5.1  11/28/95
      subroutine spline3(x,y,n,yp1,ypn,y2,u)
c***begin prologue     spline3.f
c***date written       yymmdd  
c***revision date      11/6/94      
c
c***keywords           
c***author             
c***source             @(#)spline3.f	5.1   11/28/95
c***purpose            cubic spline fit
c***description
c   given arrays x an y of length n containing a tabulated function,
c   i.e. y(i)=f(x(i)), with x1<x2<....<xn, and given values yp1 and ypn
c   for the first derivatives of the interpolating functions at points
c   1 and n, respectively, this routine returns an array y2 of length n
c   which contains the second derivative of the interpolating function 
c   at the tabulated points x(i).  if yp1 and/or ypn is greater than 10**30,
c   the routine is signalled to set the corresponding boundary condition
c   for a natural spline, with zero second derivative on that boundary.
c
c***references
c
c***routines called
c
c***end prologue       spline3.f
      implicit none
c     --- input variables -----
      integer n
      real*8 yp1,ypn
c     --- input arrays (unmodified) ---
      real*8 x(n),y(n)
c     --- input arrays (scratch) ---
      real*8 u(n)
c     --- output arrays ---
      real*8 y2(n)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i,k
      real*8 zero,half,one,two,three,six
      real*8 big,sig,p,qn,un
c
      parameter (zero=0.0d+00,half=0.5d+00,one=1.0d+00,two=2.0d+00)
      parameter (three=3.0d+00,six=6.0d+00)
      parameter (big=1.0d+30)
c
c     --- set the lower boundary condition to be 'natural' or to have
c         a prespecified first derivative.
      if(yp1.ge.big) then
         y2(1)=zero
         u(1)=zero
      else
         y2(1)=-half
         u(1)=(three/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
c
c     --- decomposition loop of the tridiagonal solution. y2 and u are
c         used for temporary storage of the decomposed factors.
      do 10 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+two
         y2(i)=(sig-one)/p
         u(i)=(six*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     $        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
   10 continue
c
c     --- set the upper boundary condition.
      if(ypn.ge.big) then
         qn=zero
         un=zero
      else
         qn=half
         un=(three/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+one)
c
c     --- the backsubstitution loop of the tridiagonal algorithm.
      do 20 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
   20 continue
c
c
      return
      end
