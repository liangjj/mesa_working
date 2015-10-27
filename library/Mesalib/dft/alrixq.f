*deck @(#)alrixq.f	5.2 2/5/95
      subroutine alrixq(p,w,n,xi)
c***begin prologue     alrixq.f
c***date written       950120  
c***revision date      2/5/95      
c
c***keywords           chebyshev,numerical integration,
c                      radial quadrature 
c***author             r.l. martin, lanl
c***source             @(#)alrixq.f	5.2   2/5/95
c***purpose            returns the points and weights associated with
c                      gaussian quadrature in the radial coordinate.
c***description
c                      this routine returns the points and weights associated
c                      with integration in the radial coordinate using the
c                      chebyshev scheme described by truetler and alrichs.
c                      as suggested by them, the chebyshev polynomials of
c                      the second kind and the mapping M4
c                      is implemented.
c                      the jacobian [r**2d(r)] is implicitly included,
c                      so the working formula for integrating a polynomial f
c                      using this routine is simply Sum f(p(i))*w(i). 
c
c                      p  ---  quadrature points (r), output.
c                      w  ---  quadrature weights, output.
c                      n  ---  number of points in the quadrature, input.
c                      xi ---  scale factor which should roughly reflect
c                              the maximum of the function to be integrated. 
c
c    
c
c***references
c                      o. treutler and r. alrichs, jcp, 103,347(1995).
c***routines called
c                      none
c***end prologue       alrixq.f
      implicit none
c
c     --- input variables ---
      integer n
      real*8 xi
c     --- input arrays (unmodified) ---
c     --- input arrays (modified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 p(n),w(n)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
      real*8 alpha,x,y,one,two
      real*8 r,dr,jacob
      parameter (alpha=0.6d+00,one=1.0d+00,two=2.0d+00)
c
      integer inp,iout
      common/io/inp,iout
c
c     --- divide the integration region into n portions using the
c         points and weights of the cheyshev polynomials of the second kind.
c         these polynomials are orthonormal on the interval (-1,1) with
c         a weight function sqrt(1-x*x). the weight function is 
c         included in the weights returned.
      call chebyq(p,w,n,2)
c
c     --- now map it onto the appropriate interval in r (0,inf).
c         use map M4, and incorporate the jacobian r**2 dr.
c         M4: r=(xi/ln2)*ln(2/(1-x)) * (1+x)**alpha
c         convenient to first put in dr term.
c         dr=[(alpha*r/(1+x)) + r/((1-x)*ln(2/(1-x))) ]dx
      y=xi/log(two)
      do 10 i=1,n
         x=p(i)
         r=(one+x)**alpha
         r=y*r*log(two/(one-x))
         dr=(alpha/(one+x)) +one/((one-x)*log(two/(one-x)))
         dr=r*dr
         jacob=r*r*dr
c
         p(i)=r
         w(i)=w(i)*jacob
   10 continue
c
c
      return
      end
