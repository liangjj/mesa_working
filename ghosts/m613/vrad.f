*deck %W% %G%
      subroutine vrad(fr,ng,nlm,mxgrd,ioff,xyzgrid,c,x,r,
     $                order,break,nbreak,ind,coef,alpha)
c***begin prologue     %M%
c***date written       940304   (yymmdd)  
c***revision date      %G%
c
c***keywords           splines, radial potential 
c***author             schneider, barry(nsf) and martin, richard(lanl)
c***source             %W% %G%
c***purpose            interpolates the spline fit for the radial potential
c                      in order to return the potential an arbritray set 
c                      of grid points 
c***description
c     
c         fr      ...  the radial function on the grid.          
c         ng      ...  the number of grid points.
c         nlm     ...  the number of spherical harmonics
c         r       ...  the radial distances on the grid.
c         order   ...  the order+1 of the spline fit.  see spline.f
c         break   ...  the break points of the spline.
c         nbreak  ...  the number of break points.
c         ind     ...  an index pointing to the nearest break point.
c         coef    ...  the spline coefficients.
c    
c
c***references
c
c***routines called
c
c***end prologue       %M%
      implicit none
c     --- input variables -----
      integer ng,nlm,mxgrd,ioff
      integer order,nbreak
      real*8 alpha
c     --- input arrays (unmodified) ---
      integer ind(ng)
      real*8 xyzgrid(mxgrd,3),c(3)
      real*8 break(nbreak+1),coef(order*nbreak,nlm)
c     --- input arrays (scratch) ---
      real*8 x(ng),r(ng)
c     --- output arrays ---
      real*8 fr(ng,nlm)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer coord,grid,lm,mflag
      integer i
      integer inp,iout
      real*8 zero,one
c
      parameter (zero=0.0d+00,one=1.0d+00)
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
c     --- find the index of the largest break point below the value r(grid).
      do 20 grid=1,ng
         call interv(break,nbreak,r(grid),ind(grid),mflag)
         if(mflag.lt.0) then
            call lnkerr('vrad: problem with interv')
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
         call ppval(r,fr(1,lm),break,coef(1,lm),ng,nbreak,order,ind,0)
c        --- multiply by 1/r to get back to the potential
         call vdiv(fr(1,lm),fr(1,lm),r,ng)
   30 continue
c
c
      return
      end
