*deck %W% %G%
      subroutine poisson(itch,d,nbf,nnp,jmat,ncoul,
     $     ndmat,nat,mxgrd,
     $     dmcut,dencut,
     $     mxgbsiz,aitch,ian,
     $     c,ex,cont,ptprim,noprim,nocont,ptcont,
     $     mxcont,nprim,ntypes,nbtype,ncont,
     $     start,nocart,nobf,maxmom,minmom,mintyp,
     $     nx,ny,nz,xyzgrid,grdwts,charge,maxl,bigl,
     $     ops,left,rnuc,amu,pwtx,rr,radii,akl,rmax,lmax,
     $     nomega,nradial,grdtyp,adjust,minesz,vorder,vlmax,
     $     vradial,vncrule)
c***begin prologue     %M%
c***date written       940304      (yymmdd)  
c***revision date      %G%      
c
c***keywords           poisson, coulomb, potential, density
c***author             martin, richard(lanl) 
c***source             %W%   %G%
c***purpose            generates the coulomb potential from the density
c***description        
c                      solves poisson equation for v, given rho
c                         (del**2) v = rho
c     
c                      this is accomplished by projecting the total
c                      density into single-center pieces and a remainder.
c                      the remainder is projected onto the atomic grids of 
c                      Becke, thereby reducing the problem to a series of 
c                      atomic poisson problems.  
c
c                      for each atom, the density is decomposed into
c                      spherical harmonic components and a radial equation
c                      is solved for the potential originating from that
c                      component. the radial equation is converted into an
c                      integral equation using the appropriate green's 
c                      function, and solved via newton-core quadrature.
c                      
c                      it should be noted that the grid used to solve the
c                      atomic problem may be  different from that used 
c                      to represent the resulting potential.
c
c***references
c
c***routines called
c
c***end prologue       %M%
      implicit none
c     --- input variables ---
      integer nbf,nnp,ncoul,ndmat,nat,mxgrd,mxgbsiz
      integer nprim,ntypes,nbtype,ncont,mxcont
      integer bigl
      integer left
      integer rmax,lmax,nomega,nradial,minesz
      integer vorder,vlmax,vradial,vncrule
      real*8 dencut,dmcut
      logical adjust
c     --- input arrays (unmodified) ---
      character*(*) ops
      character*(*) grdtyp(nat)
      integer ian(nat)
      integer ptprim(nat,ntypes),noprim(nat,ntypes),nocont(nat,ntypes)
      integer ptcont(nat,ntypes),start(nat,ntypes),nocart(0:*)
      integer nobf(ntypes),maxmom(ntypes),minmom(ntypes),mintyp(ntypes)
      integer nx(*),ny(*),nz(*)
      real*8 d(nnp,ndmat)
      real*8 c(3,nat),ex(nprim),cont(ncont)
      real*8 xyzgrid(mxgrd,3),grdwts(mxgrd,nat)
c     --- input arrays (scratch) ---
      real*8 rnuc(nat,nat),amu(nat,nat),pwtx(nat)
      real*8 rr(nat),radii(nat),akl(nat,nat)
      real*8 itch(left)
      integer aitch(*)
      integer maxl(nat)
c     --- output arrays ---
      real*8 jmat(nnp,ncoul)
      real*8 charge(nat,ndmat)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer mxgblk,ngrid,ptrad,radshls,ngb,gblksiz,grid,wts,vwts
      integer top,newgrd,iat,nugrd,nuwts,nuvwts,ztop
      integer iadtwp,wpadti
      logical prnt
      data prnt/.true./
c
      common/io/inp,iout
c      
 1060 format(8x,'grid size; atom',i3,':',16x,i6,2x,i3,' blocks')
c
c     with the changes in pm511, must now generate grid for gofish.
      mxgblk=5
      ngrid=1
      ptrad=ngrid+nat
      radshls=ptrad+rmax*nat
      ngb=radshls+nat
      gblksiz=ngb+nat
      grid=iadtwp(gblksiz+mxgblk*nat)
      wts=grid+3*mxgrd*nat
      vwts=wts+mxgrd*nat
      ztop=vwts+mxgrd*nat
      top=wpadti(ztop)
      if (top .gt. left) then
         write(iout,*) 'top,left',top,left
         call lnkerr('m511: poisson,not enough core for mkgrid')
      endif
      call mkgrid(c,ian,itch(grid),itch(wts),rmax,lmax,nomega,
     $            nradial,nat,aitch(ngrid),mxgrd,itch(vwts),rnuc,amu,
     $            pwtx,rr,adjust,radii,akl,
     $            aitch(radshls),aitch(ptrad),grdtyp)
c
c     --- redefine mxgrd to be the largest atomic grid generated 
c         this is useful since standard grids may have used less
c         than expected from the default radial and angular orders. 
      newgrd=0
      do 100 iat=1,nat
         newgrd=max(newgrd,aitch(ngrid+iat-1))
  100 continue
c
c     --- pack the points and weights arrays into (newgrd,3,nat),etc.
c         and redefine pointers
      top=wpadti(vwts)
      nugrd=iadtwp(top)
      nuwts=nugrd+3*newgrd*nat
      nuvwts=nuwts+newgrd*nat
      if(newgrd.lt.mxgrd) then
         call pakgrd(itch(grid),itch(wts),itch(vwts),mxgrd,
     $               itch(nugrd),itch(nuwts),itch(nuvwts),newgrd,nat)
         wts=grid+3*newgrd*nat
         vwts=wts+newgrd*nat
         top=wpadti(vwts+newgrd*nat)
         call vmove(itch(grid),itch(nugrd),newgrd*3*nat)
         call vmove(itch(wts),itch(nuwts),newgrd*nat)
         call vmove(itch(vwts),itch(nuvwts),newgrd*nat)
         mxgrd=newgrd
      endif
c
c     --- now that we know the size of a grid block, let's see how 
c         many grid blocks we need for each atom
c         ngb(iatom) will have number of grid blocks for atom iatom, 
c         gblksiz(block,iatom) gives number of grid points for block
c         in iatom.
      call gblk(nat,mxgbsiz,mxgblk,aitch(ngrid),
     $          aitch(ngb),aitch(gblksiz))
      if(prnt) then
         do 120 iat=1,nat
            write(iout,1060) iat,aitch(ngrid+iat-1),aitch(ngb+iat-1)
 120     continue 
      endif
c
c     --- solve poisson
      call gofish(itch(ztop),d,nbf,nnp,jmat,ncoul,
     $           ndmat,nat,mxgrd,dmcut,dencut,aitch(ngb),aitch(gblksiz),
     $           mxgblk,mxgbsiz,aitch(top),ian,c,ex,cont,ptprim,noprim,
     $           nocont,ptcont,mxcont,nprim,ntypes,nbtype,ncont,
     $           start,nocart,nobf,maxmom,minmom,mintyp,nx,ny,nz,
     $           itch(grid),itch(wts),charge,maxl,bigl,ops,
     $           left-ztop,minesz,vorder,vlmax,vradial,vncrule)
c
c
      return
      end
