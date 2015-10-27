*deck %W% %G%
      subroutine mk1grd(c,ian,grid,wts,rmax,lmax,nomega,nradial,natoms,
     $     ngrid,mxgrd,vwts,rnuc,amu,pwtx,rr,adjust,radii,usesg1,
     $     radshls,ptrad,dograd,gradwt,iatom)
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
      integer nomega,nradial,natoms,iatom
      logical adjust,usesg1,dograd
c     --- input arrays (unmodified) ---
      integer ian(natoms)
      real*8 c(3,natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 grid(mxgrd,3),wts(mxgrd),
     $     gradwt(mxgrd,3,natoms)
      real*8 radii(natoms)
      integer ngrid(natoms)
      integer radshls(natoms),ptrad(rmax,natoms)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 rnuc(*),vwts(mxgrd),amu(*),pwtx(*),rr(*)
c     --- local variables ---
      integer i
      real*8 rhomax,rhomx
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
         radii(iatom)=rhomax('slater',ian(iatom))
      endif
      if (nradial.gt. mxptsr)
     $     call lnkerr('too many radial points specified')
c
c
      i=iatom
      if(ian(i).gt.36) call lnkerr('no alpha for that yet')
      if (usesg1.and. ian(i).le.18) then
c
c           --- for now, always use 51 radial points for the 
c               standard grid.
         radshls(i)=51-1
         rhomx=rhomax('slater',ian(i))
         call sg1(grid,wts,wr,pr,rhomx,c(1,i),mxgrd,51,ian(i),
     $        ngrid,ptrad(1,i))
      else
         radshls(i)=nradial
         rhomx=rhomax('slater',ian(i))
         call atgrd(grid,wts,lmax,nomega,rmax,ngrid,
     $        pr,wr,rhomx,c(1,i),mxgrd,ptrad(1,i))
      endif
c
c     --- generate the voronoi weightings --
      call vorgrad(natoms,c,grid,wts,mxgrd,ngrid,
     $     i,vwts,rnuc,amu,pwtx,rr,adjust,radii,
     $     dograd,gradwt)
 100  continue
c
c
      return
      end
