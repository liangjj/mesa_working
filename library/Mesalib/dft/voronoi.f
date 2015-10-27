*deck @(#)voronoi.f	5.1 11/6/94
      subroutine voronoi(natoms,c,xyzgrid,wts,mxgrd,ngrid,atom,
     $                   vwts,rnuc,amu,pwtx,rr,adjust,radii)
c***begin prologue     voronoi.f
c***date written       930601  
c***revision date      11/6/94
c
c***keywords           
c***author             P. J. Hay            
c***source             @(#)voronoi.f	5.1 11/6/94
c***purpose            calculates "fuzzy" voronoi weights
c***description
c   natoms     ...     number of atoms. 
c   c          ...     atomic coordinates.
c   xyzgrid    ...     grid points at which to evaluate voronoi weights.
c   wts        ...     atomic quadrature weights.
c   mxgrd      ...     maximum number of grid points.
c   ngrid      ...     number of grid points.
c   atom       ...     the sequence number of the atom whose grid we are
c                      evaluating
c   vwts       ...     voronoi weights.
c   rnuc       ...     scratch(natoms,natoms).
c   amu        ...     scratch(natoms,natoms).
c   pwtx       ...     scratch(natoms).
c   rr         ...     scratch(natoms).
c   adjust     ...     flag which determines if the cell functions are
c                      to be scaled according to atomic radii.
c   radii      ...     atomic radii used for cell adjustment.
c
c***references
c
c                      A.D. Becke, J. Chem. Phys. 88, 2547 (1988)
c***routines called
c
c***end prologue       voronoi.f
      implicit none
c     --- input variables -----
      integer natoms,mxgrd,ngrid,atom
      logical adjust
c     --- input arrays (unmodified) ---
      real*8 c(3,natoms)
      real*8 xyzgrid(mxgrd,3)
      real*8 radii(natoms)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 vwts(mxgrd),wts(mxgrd)
c     --- output variables ---
c     --- scratch arrays ---
      real*8 rnuc(natoms,natoms),amu(natoms,natoms)
      real*8 pwtx(natoms),rr(natoms)
c     --- local variables ---
      integer igrid,k,l,kparm,kk,iatom,jatom
      real*8 xx,yy,zz,xam,yam,denom,dist
      real*8 zero,one,half,onehf
      real*8 chi,u,a
      integer iout,inp
      common /io/inp,iout
      parameter(zero=0.0d0,one=1.0d0,half=0.5d0,onehf=1.5d0)
c
c     --- check for atomic case, note early return
      if(natoms.eq.1) then
         call vfill(vwts,one,ngrid)
         return
      endif
c
c     --- calculate internuclear distances 
      do 100 iatom=1,natoms
      do 100 jatom=1,natoms
            dist=sqrt((c(1,iatom)-c(1,jatom))**2
     $     	+ (c(2,iatom)-c(2,jatom))**2  
     $    	+ (c(3,iatom)-c(3,jatom))**2)
            rnuc(iatom,jatom)=dist
 100  continue 
c
c     --- kparm is value of repeat in recursion function
      kparm=3
c
c     --- determine the hyperboloidal coordinates for this grid 
      do 210 igrid=1,ngrid
         xx=xyzgrid(igrid,1)
         yy=xyzgrid(igrid,2)
         zz=xyzgrid(igrid,3)
         do 211 k=1,natoms
            dist=(xx-c(1,k))**2 + (yy-c(2,k))**2
     $           + (zz-c(3,k))**2
            rr(k)=sqrt(dist)
 211     continue 
         do 212 k=1,natoms
            do 212 l=1,natoms
               if (k.ne.l) then
                  amu(k,l)=(rr(k)-rr(l))/rnuc(k,l)
                  if(adjust) then
                     chi=radii(k)/radii(l)
                     u=(chi-one)/(chi+one)
                     a=u/(u*u-one)
                     amu(k,l)=amu(k,l)+a*(one-amu(k,l)*amu(k,l))
                  endif
               end if
 212     continue
c
c        --- generate the cell functions
         do 220 k=1,natoms
            pwtx(k)=one
            do 230 l=1,natoms
               if(k.ne.l) then
                  xam=amu(k,l)
                  do 300 kk=1,kparm
                     yam=onehf*xam-half*xam**3
                     xam=yam
 300              continue 
                  yam=half*(one-yam)
                  pwtx(k)=pwtx(k)*yam
               end if
 230        continue
 220     continue
c
c        --- normalize cell functions
         denom=zero
         do 213 k=1,natoms
            denom=denom+pwtx(k)
 213     continue 
         vwts(igrid)=pwtx(atom)/denom
 210  continue
c
c        --- multiply voronoi weights and numerical weights
      do 240 igrid=1,ngrid
         wts(igrid)=wts(igrid)*vwts(igrid)
  240 continue
c
c
      return
      end
