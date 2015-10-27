*deck @(#)vylm.f	5.1   11/6/94
      subroutine vylm(lmax,mxgrd,ioff,xyzgrid,xyz,n,nlm,ylm,ptlm,ctheta,
     $                phi,cphi,sphi,plm,scr)
c***begin prologue     vylm.f
c***date written       940304  
c***revision date      11/6/94      
c
c***keywords           real spherical harmonics, ylm
c***author             martin, richard(lanl) 
c***source             @(#)vylm.f	5.1   11/6/94
c***purpose            generates the real spherical harmonics on a grid.
c***description
c     
c   lmax    ...  maximum angular momentum for which y(l,m) m=-l,...,+l
c               is generated.
c   mxgrd   ... leading dimension of the grid point array.
c   ioff    ... an offset into xyzgrid which tells us where the current block
c               starts
c   xyzgrid ... the grid points (mxgrd,3)
c   xyz     ... the origin.
c   n       ... number of grid points.
c   ylm     ... the real spherical harmonics (n,nlm)
c               structuring the spherical harmonics is awkward in fortan.
c               we store them as a two dimensional array. the leading index
c               refers to the point, and each column a specific l,m value.
c               access into this array is provided by a pointer array.
c   nlm     ... the number of harmonics  Sum(0,l) (2l+1) = (l+1)*(l+1)
c   ptlm    ... pointer array into ylm (0:lmax,-lmax:lmax).
c               ylm(1,l,m) can be accessed by ylm(1,ptlm(l,m))
c   ctheta  ... scratch(n)
c   phi     ... scratch(n)
c   cphi    ... scratch(n)
c   sphi    ... scratch(n)
c   plm     ... scratch(n,0:lmax)
c   scr     ... scratch(n)
c
c***references         Press,Flannery,Teukolsky,and Vetterling,
c                      Numerical Recipes, Cambridge University Press, 1986.
c
c***routines called    vlgndr,vmul
c
c***end prologue       vylm.f
      implicit none
c     --- input variables -----
      integer lmax,n,nlm,mxgrd,ioff
c     --- input arrays (unmodified) ---
      real*8 xyzgrid(mxgrd,3),xyz(3)
c     --- input arrays (scratch) ---
      real*8 ctheta(n),phi(n),cphi(n),sphi(n)
      real*8 scr(n),plm(n,0:lmax)
c     --- output arrays ---
      integer ptlm(0:lmax,-lmax:lmax)
      real*8 ylm(n,nlm)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer i,k,l,m
      integer inp,iout
      real*8 x,y,z,r
      real*8 pi,zero,one,two,three,four
      real*8 norm,mphi
      real*8 cos,sin,atan,atan2
      logical called
c
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00,three=3.0d+00)
      parameter (four=4.0d+00)
      data called/.false./
      save called,pi
c
      common/io/inp,iout
c
c     --- set up constants ---
      if(.not.called) then
         pi=four*atan(one)
         called=.true.
      endif
c
c     --- check for bad arguments
c
c     --- generate the pointer table to the ylm functions.
c         pt(l,m) tells where to find l,m in the column index.
      k=0
      do 20 l=0,lmax
         do 10 m=-l,l
            k=k+1
            ptlm(l,m)=k
   10    continue
   20 continue
c
c     --- map cartesian coordinates to spherical coordinates
c         (cos(theta) and phi) with the origin at xyz.
      do 30 i=1,n
         x=xyzgrid(i+ioff,1)-xyz(1)
         y=xyzgrid(i+ioff,2)-xyz(2)
         z=xyzgrid(i+ioff,3)-xyz(3)
         r=sqrt(x*x+y*y+z*z)
         if(r.eq.zero) then
            ctheta(i)=zero
            phi(i)=zero
         else
            ctheta(i)=z/r
            if (x.eq.zero) then
               if(y.eq.zero) then
                  phi(i)=zero
               else if(y.gt.zero) then
                  phi(i)=pi/two
               else
                  phi(i)=three*pi/two
               endif
            else 
               phi(i)=atan2(y,x)
            endif
         endif
   30 continue
c
c     --- generate the real spherical harmonics ---
      do 60 m=0,lmax
c        --- get the legendre polynomials for l=m:lmax.
c            these are the only l values consistent with this m. 
         call vlgndr(lmax,m,plm,ctheta,n,scr)
c        --- generate the phi functions
c            note the norms used for real sine,cosine functions.
         if(m.eq.0) then
            norm=one/sqrt(two*pi)
         else 
            norm=one/sqrt(pi)
         endif
         do 40 i=1,n
            mphi=float(m)*phi(i)
            cphi(i)=norm*cos(mphi)
            sphi(i)=norm*sin(mphi)
   40    continue
         do 50 l=m,lmax
            if(m.eq.0) then
               call vmul(ylm(1,ptlm(l,m)),plm(1,l),cphi,n)
            else
               call vmul(ylm(1,ptlm(l,m)),plm(1,l),cphi,n)
               call vmul(ylm(1,ptlm(l,-m)),plm(1,l),sphi,n)
            endif
   50    continue
   60 continue
c
c
      return
      end
