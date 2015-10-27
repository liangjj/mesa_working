*deck @(#)chebyq.f	5.2 2/5/95
      subroutine chebyq(p,w,n,type)
c***begin prologue     chebyq.f
c***date written       950120  
c***revision date      2/5/95      
c
c***keywords           chebyshev,numerical integration,
c                      radial quadrature 
c***author             r.l. martin, lanl
c***source             @(#)chebyq.f	5.2   2/5/95
c***purpose            returns the points and weights associated with
c                      chebyshev quadrature in the radial coordinate.
c***description
c                      this routine returns the points and weights associated
c                      with integration in the radial coordinate using the
c                      chebyshev polynomials.
c                      the weight function is implicitly included in the 
c                      weights returned, so
c                      the working formula for integrating a polynomial f(x)
c                      using this routine is simply 
c                      Int f(x)dx = Sum w(i)*f(p(i))
c
c                      p    ---  quadrature points (r), output.
c                      w    ---  quadrature weights, output.
c                      n    ---  number of points in the quadrature, input.
c                      type ---  type=1, first kind:   weight=1/sqrt(1-x**x)
c                                type=2, second kind:  weight=sqrt(1-x*x)
c
c    
c
c***references
c                      o. treutler and r. alrichs, jcp, 103,347(1995).
c***routines called
c                      none
c***end prologue       chebyq.f
      implicit none
c
c     --- input variables ---
      integer n,type
c     --- input arrays (unmodified) ---
c     --- input arrays (modified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 p(n),w(n)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
      real*8 pi,one,two,four
      real*8 x,y,arg
      data one/1.0d+00/,two/2.0d+00/,four/4.0d+00/
      save one
c
      integer inp,iout
      common/io/inp,iout
c
c
      pi=four*atan(one)
      if(type.eq.1) then
c
c        --- divide the integration region into n portions using the
c            points and weights of the cheyshev polynomials of the 
c            first kind. 
c
c            these polynomials are orthonormal on the interval (-1,1) with
c            a weight function 1/sqrt(1-x*x).
c            the general weights are associated
c            with Int f(x)/sqrt(1-x*x) dx = Sum w(i)*f(x(i)), and 
c            are w(i) = (pi/n).
c            however, we are generally dealing with an integral
c            in which the weight is not explicit.  therefore
c            we include the weight function in the weights returned; i.e.
c            define g(x)=f(x)*sqrt(1-x*x), 
c            then Int f(x)dx =Int g(x)/sqrt(1-x*x)dx = Sum w(i)*f(x(i)),
c            where w(i) = (pi/n)*sqrt(1-x*x)
c
         do 10 i=n,1,-1
            arg=float(2*i-1)/float(2*n)
            x=cos(arg*pi)
            p(i)=x
            y=pi/float(n)
            w(i)=y*sqrt(one-x*x)
   10    continue
      else if(type.eq.2) then
c
c        --- divide the integration region into n portions using the
c            points and weights of the cheyshev polynomials of the 
c            second kind.
c            these polynomials are orthonormal on the interval (-1,1) with
c            a weight function sqrt(1-x*x).  
c
c            the general weights are associated
c            with Int f(x) sqrt(1-x*x) dx = Sum w(i)*f(x(i)), and 
c            are w(i) = (pi/n+1)*sin(i*pi/n+1)*sin(i*pi/n+1).
c            however, we are generally dealing with an integral
c            in which the weight is not explicit.  therefore
c            we include the weight function in the weights returned; i.e.
c            define g(x)=f(x)/sqrt(1-x*x), 
c            then Int f(x)dx =Int g(x)*sqrt(1-x*x)dx = Sum w(i)*f(x(i)),
c            where w(i) = (pi/n+1)*sin()*sin()/sqrt(1-x*x)
c                       = (pi/n+1)*sin(i*pi/n+1)
c
         do 30 i=n,1,-1
            arg=float(i)/float(n+1)
            x=cos(arg*pi)
            p(i)=x
            y=sin(arg*pi)
            w(i)=pi*y/float(n+1)
   30    continue
      endif
c
c
      return
      end
