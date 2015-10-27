*deck @(#)splinr.f	1.1 9/8/91
c***begin prologue     splinr
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           real spline fit
c***author             schneider b (nsf)
c***source             mylib
c***purpose            spline fit
c***description        given the function and its derivative at the
c***                   two endpoints a tridiagonal set of equations
c***                   may be formulated to solve for the second
c***                   derivatives which make the interpolated second
c***                   derivatives continuous.
c
c                      x,y     = input arrays (n)
c                      yp1,ypn = input first derivatives at endpoints
c                      y2      = output second derivatives (n)
c                                (spline coefficients)
c                      u       = scratch (n)
c***references         numerical algorithms
c
c***routines called    none
c***end prologue       spline
      subroutine splinr(x,y,yp1,ypn,u,y2,n)
      implicit integer (a-z)
      real*8 x, y, yp1, ypn, u, y2, p, qn, un, sig
      real*8 test, zero, half, one, two, three, six
      dimension x(n), y(n), u(n), y2(n)
      common/io/inp,iout
      data test, zero, half, one, two, three, six / 1.d30, 0.d0, .5d0,
     1                                         1.d0, 2.d0, 3.d0, 6.d0 /
      if (yp1.gt.test) then
          y2(1)=zero
          u(1)=zero
      else
          y2(1)=-half
          u(1)=(three/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 10 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+two
         y2(i)=(sig-one)/p
         u(i)=(six*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
   10 continue
      if (ypn.gt.test) then
          qn=zero
          un=zero
      else
          qn=half
          un=(three/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
          y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+one)
      endif
      do 20 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
   20 continue
      return
      end
