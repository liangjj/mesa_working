*deck %W% %G%
      subroutine mkgrid2(c,ian,grid,wts,rmax,lmax,nomega,nradial,natoms,
     $     ngrid,mxgrd,vwts,rnuc,amu,pwtx,rr,adjust,radii,usesg1,
     $     radshls,ptrad,dograd,gradwt)
c***begin prologue     %M%
c***date written       930518  
c***revision date      %G%      
c
c***keywords           
c***author             
c***source             %W%   %G%
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       %M%
      implicit none
c     --- input variables -----
      integer rmax,lmax,mxgrd
      integer nomega,nradial,natoms
      logical adjust,usesg1,dograd
c     --- input arrays (unmodified) ---
      integer ian(natoms)
      real*8 c(3,natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 grid(mxgrd,3,natoms),wts(mxgrd,natoms),
     $     gradwt(mxgrd,3,natoms,natoms)
      real*8 radii(natoms)
      integer ngrid(natoms)
      integer radshls(natoms),ptrad(rmax,natoms)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 rnuc(*),vwts(mxgrd,natoms),amu(*),pwtx(*),rr(*)
c     --- local variables ---
      integer i
c     ---- slater's rules radii
      real*8 alpha(36)
c                h-he
      data alpha/1.00,0.59,
c                li-ne
     $           3.08,2.06, 1.55,1.23,1.04,0.89,0.77,0.68,
c                na-ar
     $           4.10,3.17, 2.59,2.17,1.89,1.66,1.47,1.33,
c                k,ca
     $           6.27,4.84, 
c                sc-zn
     $           4.59,4.38,4.20,4.01,3.82,3.69,3.53,3.40,3.27,3.16,
c                ga-kr
     $                      2.76,2.44,2.19,1.98,1.81,1.66/
      integer mxptsr
      parameter (mxptsr=250)
      real*8 pr(mxptsr),wr(mxptsr)
      integer inp,iout
      common /io/inp,iout
c
c     --- generate the atomic grids ---
      nradial=rmax-1
c     --- are we to adjust cell functions according to atomic radii.
      if(adjust) then
         do 110 i=1,natoms
            radii(i)=alpha(ian(i))
  110    continue
      endif
      if (nradial.gt. mxptsr)
     $     call lnkerr('too many radial points specified')
c
c
      do 100 i=1,natoms
         if(ian(i).gt.36) call lnkerr('no alpha for that yet')
         if (usesg1.and. ian(i).le.18) then
c
c           --- for now, always use 51 radial points for the 
c               standard grid.
            radshls(i)=51-1
            call sg1(grid(1,1,i),wts(1,i),wr,pr,alpha(ian(i)),c(1,i),
     $           mxgrd,51,ian(i),ngrid(i),ptrad(1,i))
         else
            radshls(i)=nradial
            call atgrd(grid(1,1,i),wts(1,i),lmax,nomega,rmax,ngrid(i),
     $           pr,wr,alpha(ian(i)),c(1,i),mxgrd,ptrad(1,i))
         endif
c
c        --- generate the voronoi weightings --
         call vorgrad(natoms,c,grid(1,1,i),wts(1,i),mxgrd,ngrid(i),
     $                i,vwts(1,i),rnuc,amu,pwtx,rr,adjust,radii,
     $                dograd,gradwt(1,1,1,i))
  100 continue
c
c
      return
      end
