*deck @(#)ftogrid.f	5.2 4/18/95
      subroutine ftogrid2(ng,minesz,lmax,nlm,mxgrd,ioff,xyzgrid,c,f,
     $                   fr,r,ylm,ptlm,plm,ctheta,phi,cphi,sphi,scr,
     $                   nr,rpts,ulm,y2,ind)
c***begin prologue     ftogrid.f
c***date written       940304  
c***revision date      4/18/95
c
c***keywords           
c***author             martin, richard(lanl) 
c***source             @(#)ftogrid.f	5.2 4/18/95
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
c***end prologue       ftogrid.f
      implicit none
c     --- input variables -----
      integer ng,lmax,nlm,mxgrd,ioff
      integer nr,minesz
c     --- input arrays (unmodified) ---
      real*8 xyzgrid(mxgrd,3),c(3)
      real*8 rpts(nr)
      real*8 ulm(nr,nlm),y2(nr,nlm)

c     --- input arrays (scratch) ---
      integer ind(ng)
      integer ptlm(0:lmax,0:lmax)
      real*8 fr(ng,nlm),r(ng)
      real*8 ylm(ng,nlm),plm(ng,0:lmax)
      real*8 ctheta(ng),phi(ng),cphi(ng),sphi(ng),scr(ng)
c     --- output arrays ---
      real*8 f(ng)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer lm,coord,jlo,grid,ngb,joff,gb,pt1
      integer n,highl,l,m,hnlm,oddblock
      integer inp,iout
      logical timeit
      real*8 pi,y00,vmax,vmaxl,q,vthresh
      real*8 zero,half,one,four
      real*8 timspl,timylm,tim10,timr,timind,timscrn
      real*8 dum1,dum2,dum3,dum4,dum5,dum6
c
      parameter (timeit=.false.)
      parameter (zero=0.0d0,half=0.5d0)
CCCCCC must still experiment with this flag and get it passed properly.
      data one/1.0d0/,four/4.0d0/
      data vthresh/0.0d0/
      data timspl,timylm,tim10,timr,timind,timscrn/6*0.0d0/
      save timspl,timylm,tim10,timr,timind,timscrn
      save vthresh
c
      common/io/inp,iout
c
      if(timeit) then
         call timing(dum1,dum2,dum3)
      endif
c     --- get distances on this block. use ctheta for scratch.
      call rzero(r,ng)
      do 10 coord=1,3
         call ssub(ctheta,xyzgrid(1+ioff,coord),c(coord),ng)
         call vwxy(r,r,ctheta,ctheta,+1,ng)
   10 continue
      call vsqrt(r,r,ng)
      if(timeit) then
         call timing(dum4,dum5,dum6)
         timr=timr+dum4-dum1
      endif
c
c     --- find the index of the largest break point below the value r(grid).
      if(timeit) then
         call timing(dum1,dum2,dum3)
      endif
      jlo=1
      do 20 grid=1,ng
             call hunt(rpts,nr,r(grid),jlo)
             ind(grid)=jlo
C  THIS MUST BE FIXED SO WE DONT EXTRAPOLATE
c             if(jlo.eq.0.or.jlo.eq.nr) then
c                write(iout,*) 'jlo,nr',jlo,grid,r(grid),rpts(nr)
c                call lnkerr('grid point out of spline range')
c             endif
   20 continue
      if(timeit) then
         call timing(dum4,dum5,dum6)
         timind=timind+dum4-dum1
      endif
c
c     --- put the nlm=1 (l=0,m=0) component down always.
      call splint(rpts,ulm(1,1),y2(1,1),nr,r,fr(1,1),ind,ng)
c     --- multiply by 1/r to get back to the potential, and by
c         sqrt(1/4*pi) for Y(0,0)
      pi=four*atan(one)
      y00=one/sqrt(four*pi)
      call vdiv(fr(1,1),fr(1,1),r,ng)
c      call smul(f,fr(1,1),y00,ng)
      call vwxs(f,f,fr(1,1),y00,+1,ng)
c
c     --- split the block up into smaller pieces.
      ngb=ng/minesz
      oddblock=mod(ng,minesz)
      if(oddblock.ne.0) ngb=ngb+1
      joff=0
c
      do 100 gb=1,ngb
         if(gb.ne.ngb.or.oddblock.eq.0) then
            n=minesz
         else
            n=oddblock
         endif
c
c        --- get an estimate of the magnitude of the lm components.
c            the potential should fall off as r**(-l) towards infinity
c            and so generally there should be some maximum l value
c            above which we can ignore the contribution of the lm
c            components on this grid block. this may speed up the evaluation
c            of the radial interpolations and ylm evaluation.
         if(timeit) then
            call timing(dum1,dum2,dum3)
         endif
         lm=1
         highl=0
         vmaxl=zero
         do 120 l=1,lmax
            vmax=zero
            do 115 m=-l,l
               lm=lm+1
c              --- and keep track of the largest lm we should do.
               do 113 grid=1,n
                  q=ulm(ind(grid+joff)+1,lm)+ulm(ind(grid+joff),lm)
                  vmax=max(vmax,abs(half*q))
  113          continue
  115       continue
            if(vmax.ge.vthresh) highl=max(highl,l)
            vmaxl=max(vmaxl,vmax)
c            write(iout,*) 'l,vmax',l,vmax
  120    continue
c         write(iout,*) 'highl',highl,vmaxl
         if(timeit) then
            call timing(dum4,dum5,dum6)
            timscrn=timscrn+dum4-dum1
         endif
c
c        --- do the components
         if(timeit) then
            call timing(dum1,dum2,dum3)
         endif
         lm=1
         do 140 l=1,highl
            do 135 m=-l,l
               lm=lm+1
               call splint(rpts,ulm(1,lm),y2(1,lm),nr,r(1+joff),
     $                     fr(1+joff,lm),ind(1+joff),n)
               call vdiv(fr(1+joff,lm),fr(1+joff,lm),r(1+joff),n)
  135       continue
  140    continue
         if(timeit) then
            call timing(dum4,dum5,dum6)
            timspl=timspl+dum4-dum1
         endif
c
c        --- determine ylm on the grid.
         if(timeit) then
            call timing(dum1,dum2,dum3)
         endif
         hnlm=(highl+1)*(highl+1)
         call vylm(highl,mxgrd,ioff+joff,xyzgrid,c,n,
     $             hnlm,ylm,ptlm,ctheta,phi,cphi,sphi,plm,scr)
         if(timeit) then
            call timing(dum4,dum5,dum6)
            timylm=timylm+dum4-dum1
         endif
c
c        --- generate the function
         if(timeit) then
            call timing(dum1,dum2,dum3)
         endif
c        fix this sleaze with ylm later. have to skip y00 term in frist n
c        spots.
         pt1=n
         do 99 lm=2,hnlm
            call vwxy(f(1+joff),f(1+joff),fr(1+joff,lm),
     $                ylm(1+pt1,1),+1,n)
c                     ylm(1,lm),+1,n)
            pt1=pt1+n
   99    continue
         if(timeit) then
            call timing(dum4,dum5,dum6)
            tim10=tim10+dum4-dum1
         endif
         joff=joff+n
  100 continue
c
c     --- report timings
      if(timeit) then
         write(iout,*) 'timr,timind,timscrn,timspl,timylm,tim10',
     $                  timr,timind,timscrn,timspl,timylm,tim10
      endif
c
c
      return
      end
