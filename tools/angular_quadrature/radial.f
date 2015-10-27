*deck %W% %G%
      subroutine radial(p,w,n,alpha)
      implicit none
c***begin prologue     %M%
c***date written       930501  
c***revision date      %G%      
c
c***keywords           euler-mclaurin,numerical integration,
c                      radial quadrature 
c***author             r.l. martin, lanl
c***source             %W%   %G%
c***purpose            returns the points and weights associated with
c                      gaussian quadrature in the radial coordinate.
c***description
c                      this routine returns the points and weights associated
c                      with integration in the radial coordinate using the
c                      euler-mclaurin scheme described by murray, et al.
c                      as suggested by them, the transformation with mr=2
c                      is implemented, so that the integrand and its 
c                      derivatives up to the 3m-1=5th vanish at q=0(r=0).
c                      the jacobian [r**2d(r)] is implicitly included, so
c                      the working formula for integrating a polynomial f
c                      using this routine is simply Sum f(p(i))*w(i). also
c                      note that the contributions from the function at the
c                      endpoints (i=0, i=n) are zero, so the sum runs
c                      from i=1,n-1.
c
c                      p  ---  quadrature points (r), output.
c                      w  ---  quadrature weights, output.
c                      n  ---  number of points in the quadrature, input.
c                  alpha  ---  scale factor which should roughly reflect
c                              the maximum of the function to be integrated. 
c
c    
c
c***references
c                      c.w.murray, n.c.handy, and g.j.laming, mol. phys.
c                        78,997(1993).
c***routines called
c                      none
c***end prologue       %M%
c
c     --- input variables ---
      integer n
      real*8 alpha
c     --- input arrays (unmodified) ---
c     --- input arrays (modified) ---
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 p(n),w(n)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i
      real*8 two
      parameter (two=2.0d+00)
c
      integer inp,iout
      common/io/inp,iout
c
c     --- divide the integration region into n portions --- 
c         because of the nature of the euler-mclaurin method, the integrand
c         will not contribute to the integral at the endpoints i=0 and
c         i=n, so reflect this in the points and weights returned.
      do 10 i=1,n-1
         p(i)=alpha*float(i*i)/((n-i)*(n-i))
         w(i)=p(i)*p(i)*two*alpha*float(n*i)/float((n-i)**3)
   10 continue
c
c
      return
      end
