      subroutine spline(x,y,n,yp1,ypn,scr,y2)
      implicit real*8(a-h,o-z)
      real*8 y,yp1,ypn,scr,y2,p,qn,un,sig
      dimension x(n),y(n),scr(n),y2(n)
      y2(1)=-.5
      scr(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      do 11 i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
11    scr(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1 /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*scr(i-1))/p
      qn=.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      y2(n)=(un-qn*scr(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+scr(k)
12    continue
      return
      end

