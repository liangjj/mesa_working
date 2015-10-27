*deck rk.f 
c***begin prologue     rk
c***date written       980924   (yymmdd)
c***revision date               (yymmdd)
c***keywords           odeint
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
c***end prologue       rk
c
      subroutine rk(y,dy,x,t,n,neq,nstp)
c 
      implicit integer (a-z)
      real*8 x, y, dy, h, t, dum
      dimension y(neq,nstp), dy(neq), x(nstp), t(*)
      common/io/ inp, iout
      do 10 i=1,nstp-1
         call df(x(i),y(1,i),dy,dum,idum,neq)
         h=x(i+1)-x(i)
         call rk4(y(1,i),dy,x(i),y(1,i+1),t,h,neq)
 10   continue
      return
      end


