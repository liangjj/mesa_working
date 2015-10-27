*deck @(#)mkgrid.f	5.4 11/28/95
      subroutine mkgrid(c,ian,grid,wts,rmax,lmax,nomega,nradial,natoms,
     $                  ngrid,mxgrd,vwts,rnuc,amu,pwtx,rr,adjust,radii,
     $                  akl,radshls,ptrad,grdtyp)
c***begin prologue     mkgrid.f
c***date written       930518  
c***revision date      11/28/95      
c
c***keywords           
c***author             
c***source             @(#)mkgrid.f	5.4   11/28/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       mkgrid.f
      implicit none
c     --- input variables -----
      integer rmax,lmax,mxgrd
      integer nomega,nradial,natoms
      logical adjust
c     --- input arrays (unmodified) ---
      integer ian(natoms)
      real*8 c(3,natoms)
      character*(*) grdtyp(natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 grid(mxgrd,3,natoms),wts(mxgrd,natoms)
      real*8 radii(natoms),vwts(mxgrd,natoms)
      integer ngrid(natoms)
      integer radshls(natoms),ptrad(rmax,natoms)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 rnuc(*),amu(*),pwtx(*),rr(*)
      real*8 akl(*)
c     --- local variables ---
      integer iatom
      integer inp,iout
      common /io/inp,iout
c
c     --- generate the atomic grids ---
      do 100 iatom=1,natoms
         call mkatmg(c,ian,grid(1,1,iatom),wts(1,iatom),
     $               rmax,lmax,nomega,nradial,
     $               natoms,ngrid(iatom),mxgrd,vwts(1,iatom),rnuc,amu,
     $               pwtx,rr,adjust,radii,akl,grdtyp(iatom),
     $               radshls(iatom),ptrad(1,iatom),.false.,vwts,iatom)
  100 continue
c
c
      return
      end
