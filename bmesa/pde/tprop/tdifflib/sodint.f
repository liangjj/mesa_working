*deck sodint.f 
c***begin prologue     sodint
c***date written       980924   (yymmdd)
c***revision date               (yymmdd)
c***keywords           sodint
c***                   
c***author             schneider, barry
c***source             
c***purpose            integrate a set of coupled differential equations
c***                   using the simple second order difference method
c
c***description        
c***references       
c
c***routines called   
c***end prologue       sodint
c
      subroutine sodint(y,dy,x,x1,x2,t,h,n,neq,nstp)
      real*8 y, dy, x, t, h
      real*8 x1, x2, xa, xb, fac, dum
      dimension y(neq,*), dy(*), x(*), t(*)
      common/io/ inp, iout
c     
c     to calculate the value of the y at the second point we use
c     the runge-kutta formula.  this allows us to have an entirely
c     self-starting method.  once we have two values of y we can
c     proceed with the simple second order propagation.
c
      xa=x1
      xb=xa+h
      call rk(y,dy,x,t,xa,xb,h,n,neq,2)
      fac=2.d0*h
      do 10 i=3,nstp
         x(i)=x(i-1)+h
         call df(x(i-1),y(1,i-1),dy,dum,idum,neq)
         do 20 j=1,neq
            y(j,i) = y(j,i-2) + fac*dy(j)
 20      continue
 10   continue            
         
      return
      end
