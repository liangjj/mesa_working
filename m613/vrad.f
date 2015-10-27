*deck @(#)vrad.f	5.1 11/6/94
      subroutine vrad(fr,ng,nlm,mxgrd,ioff,xyzgrid,c,x,r,
     $                nr,rpts,ulm,y2,ind)
c***begin prologue     vrad.f
c***date written       940304   (yymmdd)  
c***revision date      11/6/94
c
c***keywords           splines, radial potential 
c***author             schneider, barry(nsf) and martin, richard(lanl)
c***source             @(#)vrad.f	5.1 11/6/94
c***purpose            interpolates the spline fit for the radial potential
c                      in order to return the potential an arbritray set 
c                      of grid points 
c***description
c     
c         fr      ...  the radial function on the grid.          
c         ng      ...  the number of grid points.
c         nlm     ...  the number of spherical harmonics
c         mxgrd   ...  the maximum grid block size
c         ioff    ...  a pointer into the molecular grid array
c         xyzgrid ...  the molecular grid onto which to represent J
c         c       ...  the coordiantes of the atomic center about which the
c                      potential has been evaluated.
c         x       ...  scratch
c         r       ...  scratch
c         nr      ...  the number of radial points used in the Poisson solver.
c         rpts    ...  the radial points used in the Poisson solver.
c         ulm     ...  the function at the knots(radial points).
c         y2      ...  second derivative of interpolating function at knots. 
c         ind     ...  scratch :an index pointing to the nearest break point.
c    
c
c***references
c
c***routines called
c
c***end prologue       vrad.f
      implicit none
c     --- input variables -----
      integer ng,nlm,mxgrd,ioff
      integer nr
c     --- input arrays (unmodified) ---
      integer ind(ng)
      real*8 xyzgrid(mxgrd,3),c(3)
      real*8 rpts(nr),ulm(nr,nlm),y2(nr,nlm)
c     --- input arrays (scratch) ---
      real*8 x(ng),r(ng)
c     --- output arrays ---
      real*8 fr(ng,nlm)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer coord,grid,lm,mflag
      integer inp,iout
      real*8 zero,one
      logical debug
c
      parameter (zero=0.0d+00,one=1.0d+00)
      parameter (debug=.true.)
c 
      common/io/inp,iout
c
c     --- determine distance from the center to the grid points.
      call rzero(r,ng)
      do 10 coord=1,3
         call ssub(x,xyzgrid(1+ioff,coord),c(coord),ng)
         call vwxy(r,r,x,x,+1,ng)
   10 continue
      call vsqrt(r,r,ng)
c
c
c     --- find the index of the largest break point below the value r(grid).
      do 20 grid=1,ng
         call interv(rpts,nr-1,r(grid),ind(grid),mflag)
         if(mflag.lt.0) then
            call lnkerr('problem with spline interv')
         endif
   20 continue
c
c     --- those points which lie outside the largest break point would
c         require an extrapolation of the spline fit.  this is bad news.
c         therefore, for those points use the value of the forward 
c         integral for the evaluation.(not yet implemented).
c
c     --- get value of the radial function by interpolating the spline fit.
      do 30 lm=1,nlm
         call splint(rpts,ulm(1,lm),y2(1,lm),nr,r,fr(1,lm),ind,ng)
c        --- multiply by 1/r to get back to the potential
         call vdiv(fr(1,lm),fr(1,lm),r,ng)
   30 continue
c
c
      return
      end
