*deck %W% %G%
      subroutine ftogrid(ng,lmax,nlm,mxgrd,ioff,xyzgrid,c,f,fr,ylm,
     $                   ptlm,plm,ctheta,phi,cphi,sphi,scr,
     $                   nbreak,order,coef,break,ind,alpha)
c***begin prologue     %M%
c***date written       940304  
c***revision date      %G%
c
c***keywords           
c***author             martin, richard(lanl) 
c***source             %W% %G%
c***purpose            reconstructs the potential at an arbritrary set of
c                      points from its decomposition into spherical harmonic
c                      components and a spline fitted radial component
c                      relative to some center with coordinates c.
c***description
c                f(n) = Sum(l,m) fr(n,nlm)*ylm(n,nlm)
c   input:
c   ng      ...  the number of grid points.
c   lmax    ...  maximum l value used in spherical harmonic decomposition
c   nlm     ...  number of lm components associated with lmax.
c   ioff    ...  an offset into xyzgrid which tells us where the 
c                current grid block begins.
c   xyzgrid ...  (mxgrd,3) cartesian coordinates of the grid.
c   c            (3)    cartesian coordinates of center where expansion
c                       took place.
c   fr      ...  (ng,nlm)  radial lm component on the grid.
c   nbreak  ...  number of break points in spline fit
c   order   ...  order+1 of spline fit.
c   break   ...  spline breaks.
c   coef    ...  spline coefficients.
c
c   output:
c   f       ...  (ng)  the function to be generated.
c
c   scratch:
c   ylm     ...  (ng,nlm)  real spherical harmonics on the grid.
c   ptlm    ...   (0:lmax,0:lmax)  pointer to lm components
c   plm     ...  (ng,0:lmax)  legendre polynomials
c   ctheta  ...  (ng) 
c   phi     ...  (ng) 
c   cphi    ...  (ng) 
c   sphi    ...  (ng) 
c   scr     ...  (ng) 
c   ind     ...  (ng)
c***references
c
c***routines called
c
c***end prologue       %M%
      implicit none
c     --- input variables -----
      integer ng,lmax,nlm,mxgrd,ioff
      integer nbreak,order
      real*8 alpha
c     --- input arrays (unmodified) ---
      real*8 xyzgrid(mxgrd,3),c(3)
      real*8 break(nbreak+1)
c     --- input arrays (scratch) ---
      integer ind(ng)
      integer ptlm(0:lmax,0:lmax)
      real*8 fr(ng,nlm)
      real*8 coef(order*nbreak,nlm)
      real*8 ylm(ng,nlm),plm(ng,0:lmax)
      real*8 ctheta(ng),phi(ng),cphi(ng),sphi(ng),scr(ng)
c     --- output arrays ---
      real*8 f(ng)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer lm
      integer inp,iout
c
      common/io/inp,iout
c
c     --- determine the radial component on the grid. 
      call vrad(fr,ng,nlm,mxgrd,ioff,xyzgrid,c,ctheta,phi,
     $          order,break,nbreak,ind,coef,alpha)
c
c     --- determine ylm on the grid.
      call vylm(lmax,mxgrd,ioff,xyzgrid,c,ng,
     $          nlm,ylm,ptlm,ctheta,phi,cphi,sphi,plm,scr)
c
c     --- generate the function
      do 10 lm=1,nlm
         call vwxy(f,f,fr(1,lm),ylm(1,lm),+1,ng)
   10 continue
c
c
      return
      end
