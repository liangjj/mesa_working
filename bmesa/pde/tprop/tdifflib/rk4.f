*deck rk4.f 
c***begin prologue     rk4
c***date written       980924   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   
c***author             numerical recipe
c***source             
c***purpose            integrate a set of coupled differential equations
c***                   using a fourth order runge-kutta method.  
c
c***description        method is described in numerical recipes
c***references       
c
c***routines called   
c***end prologue       rk4
c
      subroutine rk4(y,dy,x,yout,t,h,n)
c 
      implicit integer (a-z)
      real*8 x, y, dy, yout, h, h6, hh, xh, t, dum
      dimension y(n), dy(n), yout(n), t(n,3)
      common/io/ inp, iout
      hh=h*.5d0
      h6=h/6.d0
      xh=x+hh
      do 10 i=1,n
         t(i,1) = y(i) + hh*dy(i)
 10   continue
      call df(xh,t(1,1),t(1,2),dum,idum,n)
      do 20 i=1,n
          t(i,1) = y(i) + hh*t(i,2)
 20   continue
      call df(xh,t(1,1),t(1,3),dum,idum,n)
      do 30 i=1,n
         t(i,1) = y(i) + h*t(i,3)
         t(i,3) = t(i,2) + t(i,3)
 30   continue
      call df(x+h,t(1,1),t(1,2),dum,idum,n)
      do 40 i=1,n
         yout(i) = y(i) + h6*( dy(i) + t(i,2) + 2.d0*t(i,3) )                             
 40   continue                   
      return
      end


