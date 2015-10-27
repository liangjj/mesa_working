*deck @(#)getgrid.f	5.1  4/18/95
      subroutine getgrid(itch,aitch,left,nat,ian,
     $     c,ops,rmax,lmax,grdtyp,adjust,
     $     mxgrd,ngrid,grid,wts,top)
c***begin prologue     getgrid.f
c***date written       yymmdd  
c***revision date      4/18/95      
c
c***keywords           
c***author             
c***source             @(#)getgrid.f	5.1   4/18/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       getgrid.f
c     --- input variables ---
      integer left,nat,rmax,lmax
      integer top
      logical adjust
c     --- input arrays (unmodified) ---
      character*(*) ops
      character*(*) grdtyp(nat)
      integer ian(nat)
      real*8 c(3,nat)
c     --- input arrays (scratch) ---
      real*8 itch(left)
      integer aitch(*)
c     --- output arrays ---
c     --- output variables ---
      integer ngrid,grid,wts
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer nomega,nradial
      integer ptrad,radshls,akl
      integer nugrd,nuwts
      integer iadtwp,wpadti
      logical prnt
      data prnt/.true./
c
      common/io/inp,iout
c      
 1060 format(8x,'grid size; atom',i3,':',16x,i6,2x,i3,' blocks')
c
c     --- generate the integration grids. 
c         allocate some scratch memory. 
      ngrid=top
      grid=iadtwp(ngrid+nat)
      wts=grid+3*mxgrd*nat
      ptrad=wpadti(wts+mxgrd*nat)
      vwts=iadtwp(ptrad+rmax)
      rnuc=vwts+mxgrd
      amu=rnuc+nat*nat
      pwtx=amu+nat*nat
      rr=pwtx+nat
      radii=rr+nat
      akl=radii+nat
      radshls=akl+nat*nat
      top=wpadti(radshls+nat)
      if (top .gt. left) then
         write(iout,*) 'top,left',top,left
         call lnkerr('m511: poisson,not enough core for mkgrid')
      endif
      call mkgrid(c,ian,itch(grid),itch(wts),rmax,lmax,nomega,
     $            nradial,nat,aitch(ngrid),mxgrd,itch(vwts),
     $            itch(rnuc),itch(amu),itch(pwtx),itch(rr),
     $            adjust,itch(radii),itch(akl),
     $            aitch(radshls),aitch(ptrad),grdtyp)
c
c     --- now determine the largest grid size actually generated.
c         redefine mxgrd to be the largest atomic grid generated 
c         this is useful since standard grids may have used less
c         than expected from the default radial and angular orders. 
      biggrd=0
      do 10 i=1,nat
         biggrd=max(biggrd,aitch(ngrid+i-1))
   10 continue
c
c     --- pack the points and weights arrays into (biggrd,3,nat),etc.
c         and redefine pointers
      nugrd=iadtwp(top)
      nuwts=nugrd+3*biggrd*nat
      if(biggrd.lt.mxgrd) then
         call pakgrd(itch(grid),itch(wts),mxgrd,itch(nugrd),itch(nuwts),
     $               biggrd,nat)
         wts=grid+3*biggrd*nat
         top=wpadti(wts+biggrd*nat)
         call vmove(itch(grid),itch(nugrd),biggrd*3*nat)
         call vmove(itch(wts),itch(nuwts),biggrd*nat)
         mxgrd=biggrd
      endif
c
c
      return
      end
