*deck secder.f
c***begin prologue     secder
c***date written       930623   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           second derivative of ricatti-bessel function
c***author             schneider, barry (nsf)
c***source             m6202
c***                   
c***                   
c
c***references         
c
c***routines called    
c***end prologue       secder
c
      subroutine secder(f,ddf,x,l,n)
      implicit integer (a-z)
      real*8 f, ddf, x
      dimension f(n), ddf(n), x(n)
c
      do 10 i=1,n
         ddf(i)=( l*(l+1) - x(i)*x(i) )*f(i)/(x(i)*x(i))
   10 continue     
      return
      end
