*deck @(#)getgrid.f	1.2  11/28/95
      subroutine getgrid(itch,aitch,left,nat,ian,
     $                   c,ops,rmax,lmax,grdtyp,adjust,
     $                   mxgrd,ngrid,grid,wts,vwts,ptrad,maxr,
     $                   top)
c***begin prologue     getgrid.f
c***date written       yymmdd  
c***revision date      11/28/95
c
c***keywords           
c***author             
c***source             @(#)getgrid.f	1.2   11/28/95
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
      implicit none
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
      integer ngrid,grid,wts,vwts,ptrad,mxgrd,maxr
c     --- scratch arrays ---
c     --- local variables ---
      integer i
      integer inp,iout
      integer nomega,nradial,angsiz
      integer radshls,akl
      integer biggrd,nugrd,nuwts
      integer iadtwp,wpadti
      integer rnuc,amu,pwtx,rr,radii,nuvwts,nuptrad
      integer atom
      integer stderr,ierr
      logical prnt
      data prnt/.true./
c
      common/io/inp,iout
c      
 1060 format(8x,'grid size; atom',i3,':',16x,i6,2x,i3,' blocks')
      ierr=stderr()
c
c     --- generate the integration grids. 
c         allocate some scratch memory. 
      ngrid=top
      grid=iadtwp(ngrid+nat)
      wts=grid+3*mxgrd*nat
      vwts=wts+mxgrd*nat
      ptrad=wpadti(vwts+mxgrd*nat)
      rnuc=iadtwp(ptrad+rmax*nat)
      amu=rnuc+nat*nat
      pwtx=amu+nat*nat
      rr=pwtx+nat
      radii=rr+nat
      akl=radii+nat
      radshls=wpadti(akl+nat*nat)
      top=radshls+nat
      if (top .gt. left) then
         write(iout,*) 'top,left',top,left
         call plnkerr('m515: getgrid,not enough core for mkgrid',2000)
      endif
      nomega=angsiz(lmax)
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
c     --- and the largest radial block.
      maxr=0
      do 12 i=1,nat
         maxr=max(maxr,aitch(radshls+i-1))
   12 continue
c
c     --- pack the points and weights arrays into (biggrd,3,nat),etc.
c         and redefine pointers
      nugrd=iadtwp(top)
      nuwts=nugrd+3*biggrd*nat
      nuvwts=nuwts+biggrd*nat
      nuptrad=wpadti(nuvwts+biggrd*nat)
      top=nuptrad+maxr*nat
      if(top.gt.left) then
         write(iout,*) 'top,left',top,left
         call plnkerr('need more core for getgrid',2001)
      endif
      if(biggrd.lt.mxgrd) then
         call pakgrd(itch(grid),itch(wts),itch(vwts),mxgrd,
     $               itch(nugrd),itch(nuwts),itch(nuvwts),biggrd,nat)
         call icopy(aitch(ptrad),aitch(nuptrad),rmax*nat)
         wts=grid+3*biggrd*nat
         vwts=wts+biggrd*nat
         ptrad=wpadti(vwts+biggrd*nat)
         top=ptrad+maxr*nat
         call vmove(itch(grid),itch(nugrd),biggrd*3*nat)
         call vmove(itch(wts),itch(nuwts),biggrd*nat)
         call vmove(itch(vwts),itch(nuvwts),biggrd*nat)
         do 15 atom=1,nat
            do 14 i=1,maxr
               aitch(ptrad+i-1+(atom-1)*maxr)=
     $            aitch(nuptrad+i-1+(atom-1)*rmax)
   14       continue
   15    continue
         mxgrd=biggrd
      endif
c
c
      return
      end
