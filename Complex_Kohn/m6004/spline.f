*deck spline
c***begin prologue     spline
c***date written       890404   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           spline, link 6004, orbital decomposition
c***author             rescigno, t. n.(llnl)
c***source             spline
c***purpose            spline fit
c***references       
c
c***routines called    none
c***end prologue       spline
      subroutine spline(x,y,n,yp1,ypn,scr,y2)
      implicit integer (a-z)
      real *8 x
      complex*16 y, yp1, ypn, scr, y2, p, qn, un, sig
      dimension x(n), y(n), scr(n), y2(n)
      y2(1)=-.5d0
      scr(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      do 10 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.d0
         y2(i)=(sig-1.d0)/p
         scr(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*scr(i-1))/p
   10 continue
      qn=.5d0
      un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      y2(n)=(un-qn*scr(n-1))/(qn*y2(n-1)+1.)
      do 20 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+scr(k)
   20 continue
      return
      end
